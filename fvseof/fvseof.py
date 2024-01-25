import cobra
from cobra import Reaction
from cobra.flux_analysis.loopless import loopless_solution
from cobra.flux_analysis import flux_variability_analysis
import numpy as np
import pandas as pd

class FVSEOF():
    def __init__(self, model, biomass_reaction_id, target_metabolite_id, essential_reaction_threshold = 0.5):
        """
        Initialize the FVSEOF class.

        Parameters
        ----------
        model : cobra.Model
            The model to perform FVSEOF on.
        biomass_reaction_id : str
            The id of the biomass reaction in the model.
        target_metabolite_id : str
            The id of the metabolite to perform FVSEOF on.
        essential_reaction_threshold : float, optional
            The threshold fraction of max biomass production below which knocked out reactions are considered essential. The default is 0.5.
        """
        self.model = model
        self.biomass_reaction_id = biomass_reaction_id
        self.target_metabolite_id = target_metabolite_id
        assert self.biomass_reaction_id in [r.id for r in self.model.reactions], "Biomass reaction not in model."
        assert self.target_metabolite_id in [m.id for m in self.model.metabolites], "Target metabolite not in model."

        self.product_sink_reaction_id = self.add_product_sink_reaction(self.target_metabolite_id)
        self.product_max_theoretical_yield = self.calc_product_maximal_theoretical_yield()
        self.maximal_biomass_growth = self.calc_maximal_biomass_growth()

        self.check_essential_reaction_threshold = essential_reaction_threshold

    def add_product_sink_reaction(self, target_metabolite_id):
        """
        Add a sink reaction for the target metabolite to the model.

        Parameters
        ----------
        target_metabolite_id : str
            The id of the target metabolite.

        Returns
        -------
        product_sink_reaction : cobra.Reaction
            The product sink reaction.

        """
        reaction_id = target_metabolite_id+"_fvseof_sink"
        product_sink_reaction = Reaction(target_metabolite_id+"_fvseof_sink")
        product_sink_reaction.add_metabolites({
            self.model.metabolites.get_by_id(target_metabolite_id): -1
        })
        product_sink_reaction.bounds = (0.0, 1000.0)
        self.model.add_reactions([product_sink_reaction])

        return reaction_id
    
    def calc_product_maximal_theoretical_yield(self):
        """
        Calculate the maximal theoretical yield of the target metabolite.

        Parameters
        ----------
        product_sink_reaction : cobra.Reaction
            The product sink reaction.

        Returns
        -------
        max_theoretical_yield : float
            The maximal theoretical yield of the target metabolite.

        """
        with self.model as model:
            model.objective = self.product_sink_reaction_id
            sol = loopless_solution(model)
            return sol.fluxes[self.product_sink_reaction_id]
        
    def calc_maximal_biomass_growth(self):
        """
        Calculate the maximal growth when product formation is not constrained. This is required for the calc_essential_reaction function.
        
        Returns
        -------
        max_growth : float
            The maximal growth.
        """
        with self.model as model:
            model.objective = self.biomass_reaction_id
            sol = loopless_solution(model)
            return sol.fluxes[self.biomass_reaction_id]
        
    def check_essential_reaction(self, reaction_id):
        """
        Check if a reaction is essential.
        
        Parameters
        ---------- 
        reaction_id : str
            The id of the reaction to check.
            
        Returns
        -------
        essential : bool
            Whether the reaction is essential or not.
        """
        with self.model as model:
            model.objective = self.biomass_reaction_id
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.lower_bound = 0.0
            reaction.upper_bound = 0.0
            
            try:
                sol = loopless_solution(model)
                if sol.fluxes[self.biomass_reaction_id] / self.maximal_biomass_growth < self.check_essential_reaction_threshold:
                    return True
                return False
            except:
                return True

    def run(self, n_steps = 10, check_essentiality=False, fva = True, fva_n_processes = 1) -> pd.DataFrame:
        """
        Run FVSEOF on the model. 2 points for lower bound on target production are used to find fluxes that increase or decrease when the target production is increased.
        
        Parameters
        ----------
        n_steps : int, optional
            The number of steps to use for the gradual increase on the lower bound on target production. The default is 10.
        check_essentiality : bool, optional
            Whether to check if target reactions are essential. The default is True.
        fva : bool, optional
            Whether to perform flux variability analysis. If False, FBA is used, which turns this function into the FSEOF algorithm instead of the FVSEOF algorithm. Default is True.
        fva_n_processes : int, optional
            The number of processes to use for flux variability analysis. The default is 1.
            
        Returns
        -------
        targets: pd.DataFrame
            A dataframe with the targets for upregulation and downregulation.
        """
        steps_lower_bounds = [n/n_steps*self.product_max_theoretical_yield for n in range(0, n_steps)]
        method = "FVA" if fva else "FBA"

        # Perform fva / fba and retain (mean) fluxes for each step
        per_step_fluxes = {r_id: [] for r_id in [r.id for r in self.model.reactions]}
        for i, step_lower_bound in enumerate(steps_lower_bounds):
            print("\rPerforming {} for step {}/{}...".format(method, i+1, len(steps_lower_bounds)), end="")
            with self.model as model:
                model.objective = self.biomass_reaction_id
                model.reactions.get_by_id(self.product_sink_reaction_id).lower_bound = step_lower_bound

                if fva:
                    fva_df = flux_variability_analysis(model, fraction_of_optimum=0.95, processes=fva_n_processes)
                    for r_id in [r.id for r in self.model.reactions]:
                        per_step_fluxes[r_id].append((fva_df.loc[r_id, "maximum"] + fva_df.loc[r_id, "minimum"]) / 2)
                else:
                    fba_sol = loopless_solution(model)
                    for r_id in [r.id for r in self.model.reactions]:
                        per_step_fluxes[r_id].append(fba_sol.fluxes[r_id])
        print("\nDone.")

        # Assign target types based on min and max fluxes
        target_types = {}
        for r_id in [r.id for r in self.model.reactions]:
            if not per_step_fluxes[r_id][0] == 0 and not per_step_fluxes[r_id][-1] == 0:
                if per_step_fluxes[r_id][0] * per_step_fluxes[r_id][-1] >= 0:
                    if abs(per_step_fluxes[r_id][-1]) > abs(per_step_fluxes[r_id][0]):
                        target_types[r_id] = "Up"
                    else:
                        target_types[r_id] = "Down"
                elif per_step_fluxes[r_id][0] * per_step_fluxes[r_id][-1] < 0:
                    target_types[r_id] = "Reverse"

        # Find slope of linear regressions
        slopes = {}
        for r_id in [r.id for r in self.model.reactions]:
            slopes[r_id] = np.polyfit(steps_lower_bounds, per_step_fluxes[r_id], 1)[0]

        # Find essentialities
        essentialities = {}
        if check_essentiality:
            print("\nCalculating essentialities...")
            for i, r_id in enumerate([r.id for r in self.model.reactions]):
                print("\rChecking essentiality for reaction " + str(i+1) + "/" + str(len([r.id for r in self.model.reactions])) + "...", end="")
                if r_id in target_types:
                    essentialities[r_id] = self.check_essential_reaction(r_id) if check_essentiality else None
                    
        # Combine into big df
        df_data = {}
        for r_id in [r.id for r in self.model.reactions]:
            if r_id in target_types:
                df_data[r_id] = {
                    "target_type": target_types[r_id],
                    "slope": slopes[r_id],
                    "gene_reaction_rule": self.model.reactions.get_by_id(r_id).gene_reaction_rule
                }

                if check_essentiality:
                    df_data[r_id]["essentiality"] = essentialities[r_id]

                for i, step_lower_bound in enumerate(steps_lower_bounds):
                    df_data[r_id]["step_" + str(i)] = per_step_fluxes[r_id][i]
            
        df = pd.DataFrame.from_dict(df_data, orient="index")
        df.index.name = "reaction_id"
        df = df.sort_values(by=['target_type', 'slope'], key=lambda x: x.abs() if x.name == 'slope' else x, ascending=[False, False])
                       
        return df