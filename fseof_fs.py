import cobra
from cobra import Reaction
from cobra.flux_analysis.loopless import loopless_solution
import numpy as np
import tqdm
import pandas as pd
import logging

class FSEOF_FS():
    def __init__(self, model, biomass_reaction_id, target_metabolite_id, essential_reaction_threshold = 0.5):
        """
        Initialize the FSEOF_FS class.

        Parameters
        ----------
        model : cobra.Model
            The model to perform FSEOF_FS on.
        biomass_reaction_id : str
            The id of the biomass reaction in the model.
        target_metabolite_id : str
            The id of the metabolite to perform FSEOF_FS on.
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

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)

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
        reaction_id = target_metabolite_id+"_sink"
        product_sink_reaction = Reaction(target_metabolite_id+"_sink")
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

    def get_average_fluxes_from_n_solutions(self, n):
        """
        Get the average fluxes from n solutions.
        
        Parameters
        ----------
        n : int
            The number of samples for flux sampling
            
        Returns
        -------
        flux_list : dict
            A dictionary with the average fluxes for each reaction.
        """
        flux_list = {}
        for _ in tqdm.tqdm(range(n)):
            sol = loopless_solution(self.model)
            for r in self.model.reactions:
                if r.id in flux_list:
                    flux_list[r.id].append(sol.fluxes[r.id])
                else:
                    flux_list[r.id] = [sol.fluxes[r.id]]
    
        for r, fs in flux_list.items():
            flux_list[r] = sum(fs)/len(fs)
    
        return flux_list

    def run(self, fraction_low=0.1, fraction_high=0.2, check_essentiality=True, n=100) -> pd.DataFrame:
        """
        Run FSEOF_FS on the model. 2 points for lower bound on target production are used to find fluxes that increase or decrease when the target production is increased.
        
        Parameters
        ----------
        fraction_low : float, optional
            The fraction of maximal theoretical yield to use for the first lower bound. The default is 0.1.
        fraction_high : float, optional
            The fraction of maximal theoretical yield to use for the second lower bound. The default is 0.2.
        check_essentiality : bool, optional
            Whether to check if target reactions are essential. The default is True.
        n : int
            The number of samples for flux sampling
            
        Returns
        -------
        targets: pd.DataFrame
            A dataframe with the targets for upregulation and downregulation.
        """
        # Fluxes at low lower bound on maximal theoretical yield production
        with self.model as model:
            model.objective = self.biomass_reaction_id
            self.model.reactions.get_by_id(self.product_sink_reaction_id).lower_bound = self.product_max_theoretical_yield * fraction_low
            self.logger.info("Getting average fluxes from " + str(n) + " samples for lower point target production...")
            fluxes_low = self.get_average_fluxes_from_n_solutions(n)

        # Fluxes at high lower bound on maximal theoretical yield production
        with self.model as model:
            model.objective = self.biomass_reaction_id
            self.model.reactions.get_by_id(self.product_sink_reaction_id).lower_bound = self.product_max_theoretical_yield * fraction_high
            self.logger.info("Getting average fluxes from " + str(n) + " samples for higher point target production...")
            fluxes_high = self.get_average_fluxes_from_n_solutions(n)

        self.logger.info("Running FSEOF_FS (check_essentiality = " + str(check_essentiality) + ")...")
        df_data = []
        for r_id in tqdm.tqdm([r.id for r in self.model.reactions]):
            fluxes_low_r = fluxes_low[r_id]
            fluxes_high_r = fluxes_high[r_id]

            mean_flux_low_r = np.mean(fluxes_low_r)
            mean_flux_high_r = np.mean(fluxes_high_r)

            # If same sign
            if mean_flux_low_r * mean_flux_high_r > 0:
                # If increased
                if abs(mean_flux_high_r) > abs(mean_flux_low_r):
                    target_type = "Up"
                # If decreased
                else:
                    target_type = "Down"

                reaction_data = {
                    "target_type": target_type,
                    "reaction_id": r_id,
                    "reaction_name": self.model.reactions.get_by_id(r_id).name,
                    "mean_flux_low": mean_flux_low_r,
                    "mean_flux_high": mean_flux_high_r
                }

                if check_essentiality:
                    essential = self.check_essential_reaction(r_id)
                    reaction_data["essential"] = essential

                df_data.append(reaction_data)

        df = pd.DataFrame(df_data)
        
        # Sort
        df["abs_diff_mean_fluxes"] = df.apply(lambda x: abs(x["mean_flux_low"] - x["mean_flux_high"]), axis=1)
        df = df.sort_values(by=["target_type", "abs_diff_mean_fluxes"], ascending=[False, False])

        return df