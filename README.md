# FSEOF-FS
A python implementation of the Flux Sampling based on Enforced Objective Flux (FSEOF) algorithm for metabolic engineering (by Choi et al. 2010) extended with Flux Sampling. In FSEOF, Flux Balance Analysis (FBA) is performed with biomass as objective function, under the constraint of a step-wise increasing lower bound on production of the target metabolite. Reaction fluxes that are positively correlated to increasing lower-bound on target metabolite production are identified as targets for upregulation in metabolic engineering, while reaction fluxes that are negatively correlated to increasing lower-bound on target metabolite production are identified as targets for downregulation.

In FVSEOF (Park et al. 2012), Flux Variability Analysis was proposed to address the non-uniqueness of FBA by evaluating the possible range of fluxes in the optimal solution space. In this implementation, Flux Sampling was chosen as alternative to FVA, as it also yields approximate flux ranges matching the optimal solution, but requires much less computational time. 



## References
 
* FSEOF: Choi, H. S., Lee, S. Y., Kim, T. Y., & Woo, H. M. (2010). In Silico Identification of Gene Amplification Targets for Improvement of Lycopene Production. In Applied and Environmental Microbiology (Vol. 76, Issue 10, pp. 3097–3105). American Society for Microbiology. https://doi.org/10.1128/aem.00115-10
* FVSEOF: Park, J. M., Park, H. M., Kim, W. J., Kim, H. U., Kim, T. Y., & Lee, S. Y. (2012). Flux variability scanning based on enforced objective flux for identifying gene amplification targets. In BMC Systems Biology (Vol. 6, Issue 1). Springer Science and Business Media LLC. https://doi.org/10.1186/1752-0509-6-106
