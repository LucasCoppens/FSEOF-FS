# FVSEOF
Flux Variability Sampling based on Enforced Objective Flux (FVSEOF, Park et al. 2012) is an algorithm that uses Genome-scale metabolic models (GEMs) to find gene targets for upregulation and downregulation that increase metabolic flux towards a desired metabolic product.

This package provides a python implementation of the FVSEOF algorithm. In FVSEOF, Flux Variability Analysis (FVA) is performed with biomass as objective function, under the constraint of a step-wise increasing lower bound on production of the target metabolite. Reaction fluxes that are positively correlated to increasing lower-bound on target metabolite production are identified as targets for upregulation in metabolic engineering, while reaction fluxes that are negatively correlated to increasing lower-bound on target metabolite production are identified as targets for downregulation.

When running the algorithm it is also possible to use Flux Balance Analysis (FBA) instead of FVA. This effectively turns the algorithm into the earlier version FSEOF by Choi et al. (2010). This option is a lot faster but may yield different results. 

This Python implementation of FVSEOF and FSEOF is based on the work published by [Park et al. (2012)](https://doi.org/10.1186/1752-0509-6-106) and [Choi et al. (2010)](https://doi.org/10.1128/aem.00115-10). All credit for the design and concept of the algorithm goes to the original authors.

## Installation

```bash
# 1. Clone the repository
git clone git@github.com:LucasCoppens/fvseof.git

# 2. Change directory to the cloned repository
cd fvseof

# 3. install using python
python3 -m pip install .
```

## Usage

### Basic

```python
# Import the FSEOF_FS class and cobra
from fvseof.fvseof import FVSEOF
import cobra

# Specify model sbml filename, biomass reaction ID and target metabolite ID
model_filename = "path/to/GEM"
model_biomass_reaction_id = "bio1" 
model_target_metabolite_id = "cpd00029_c"

# Use cobra to read sbml file
model = cobra.io.read_sbml_model(model_filename)
model.solver = "glpk"

# Run fseof_fs
fvseof = FVSEOF(model, model_biomass_reaction_id, model_target_metabolite_id)
df = fvseof.run()
```

### Varying parameters

```python
fvseof = FVSEOF(model, model_biomass_reaction_id, model_target_metabolite_id)

# Using an amount of steps other than 10
df = fvseof.run(n_steps = 15)

# Use parallel processes to speed up FVA
df = fvseof.run(fva_n_processes = 8)

# Return the df with an "essentiality" column which has information on whether a reaction is essential for biomass production in the model
df = fvseof.run(check_essentiality=True)

# Running FSEOF instead of FVSEOF
df = fvseof.run(fva=False)
```

## References
 
* FVSEOF: Park, J. M., Park, H. M., Kim, W. J., Kim, H. U., Kim, T. Y., & Lee, S. Y. (2012). Flux variability scanning based on enforced objective flux for identifying gene amplification targets. In BMC Systems Biology (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1752-0509-6-106
* FSEOF: Choi, H. S., Lee, S. Y., Kim, T. Y., & Woo, H. M. (2010). In Silico Identification of Gene Amplification Targets for Improvement of Lycopene Production. In Applied and Environmental Microbiology (Vol. 76, Issue 10, pp. 3097–3105). American Society for Microbiology. https://doi.org/10.1128/aem.00115-10
