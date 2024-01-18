# FSEOF_FS
Flux Sampling based on Enforced Objective Flux (FSEOF, Choi et al. 2010) is an algorithm that uses Genome-scale metabolic models (GEMs) to find gene targets for upregulation and downregulation that increase metabolic flux towards a desired metabolic product.

This package provides a python implementation of the FSEOF algorithm, extended with Flux Sampling (FSEOF_FS). In FSEOF, Flux Balance Analysis (FBA) is performed with biomass as objective function, under the constraint of a step-wise increasing lower bound on production of the target metabolite. Reaction fluxes that are positively correlated to increasing lower-bound on target metabolite production are identified as targets for upregulation in metabolic engineering, while reaction fluxes that are negatively correlated to increasing lower-bound on target metabolite production are identified as targets for downregulation.

In FVSEOF (Park et al. 2012), Flux Variability Analysis was proposed to address the non-uniqueness of FBA by evaluating the possible range of fluxes in the optimal solution space. In this implementation, Flux Sampling was chosen as alternative to FVA, as it also yields approximate flux ranges matching the optimal solution, but requires much less computational time. 

In FSEOF_FS, two sets of solutions are calculated by performing flux sampling on the model under two different lower bound constraints on the production of the desired product. These lower bounds are set as fractions of the maximal theoretical yield of the product, which is first calculated by the algorithm. The default constraints on product formation lower bound are 10% (0.1) and 90% (0.9), but they can be set by the user. The biomass reaction is set as the objective. After flux sampling, reactions with mean fluxes of the same sign between the two sets are considered. If mean flux of a reaction increases with increasing lower bound on product formation, it is identified as a target for upregulation. If mean flux of a reaction decreases with increasing lower bound on product formation, it is identified as a target for downregulation.

## Installation

```bash
# 1. Clone the repository
git clone git@github.com:LucasCoppens/fseof_fs.git

# 2. Change directory to the cloned repository
cd fseof_fs

# 3. install using python
python3 -m pip install .
```

## Usage

### Basic

```python
# Import the FSEOF_FS class and cobra
from fseof_fs.fseof_fs import FSEOF_FS
import cobra

# Specify model sbml filename, biomass reaction ID and target metabolite ID
model_filename = "path/to/GEM"
model_biomass_reaction_id = "bio1" 
model_target_metabolite_id = "cpd00029_c"

# Use cobra to read sbml file
model = cobra.io.read_sbml_model(model_filename)
model.solver = "glpk"

# Run fseof_fs
fseof_fs = FSEOF_FS(model, model_biomass_reaction_id, model_target_metabolite_id)
df = fseof_fs.run(check_essentiality=False, n=100, fraction_low=0.1, fraction_high=0.9)
```


## References
 
* FSEOF: Choi, H. S., Lee, S. Y., Kim, T. Y., & Woo, H. M. (2010). In Silico Identification of Gene Amplification Targets for Improvement of Lycopene Production. In Applied and Environmental Microbiology (Vol. 76, Issue 10, pp. 3097–3105). American Society for Microbiology. https://doi.org/10.1128/aem.00115-10
* FVSEOF: Park, J. M., Park, H. M., Kim, W. J., Kim, H. U., Kim, T. Y., & Lee, S. Y. (2012). Flux variability scanning based on enforced objective flux for identifying gene amplification targets. In BMC Systems Biology (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1752-0509-6-106
