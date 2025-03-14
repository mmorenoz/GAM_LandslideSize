#  Modeling the area of co-seismic landslides via data-driven models: The Kaikōura example
This is the R code to reproduce the analysis in "Modeling the area of co-seismic landslides via data-driven models: The Kaikōura example" by 
Mateo Moreno <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-9530-3076)</sup>, 
Stefan Steger <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0886-5191)</sup>, 
Hakan Tanyas <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-0609-2140)</sup> and 
Luigi Lombardo <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-4348-7288)</sup>.

> Moreno, M., Steger, S., Tanyas, H., & Lombardo, L. (2023). Modeling the area of co-seismic landslides via data-driven models: The Kaikōura example. Engineering Geology, 320, 107121. https://doi.org/10.1016/j.enggeo.2023.107121

The provided code runs the landslide area modeling using GAMS, its respective validation via random cross-validation and spatial cross-validation, as well as the visualization of the results.


## Abstract
The last three decades have witnessed substantial developments in data-driven models for landslide prediction. However, this improvement has been mostly devoted to models that estimate locations where landslides may occur, i.e., landslide susceptibility. Although susceptibility is crucial when assessing landslide hazard, another equally important piece of information is the potential landslide area once landslides initiate on a given slope. This manuscript addresses this gap in the literature by using a Generalized Additive Model whose target variable is the topographically-corrected landslide areal extent at the slope unit level. In our case, the underlying assumption is that the variability of landslide areas across the geographic space follows a Log-Normal probability distribution. We test this framework on co-seismic landslides triggered by the Kaikōura earthquake (7.8 Mw on November 13th 2016). The model performance was evaluated using random and spatial cross-validation. Additionally, we simulated the expected landslide area over slopes in the entire study area, including those that did not experience slope failures in the past. The performance scores revealed a moderate to strong correlation between the observed and predicted landslide areas. Moreover, the simulations show coherent patterns, suggesting that it is worth extending the landslide area prediction further. We share data and codes in a GitHub repository to promote the repeatability and reproducibility of this research.

## Repo structure
The general structure is as follows:
- `dat`: data sets
- `dev`: development (scripts)

