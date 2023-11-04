# Global disparities in COVID-19 vaccine coverage associated with trajectories of SARS-CoV-2 adaptation

# Calculation of the ratio of nonsynonymous to synonymous divergence （dN/dS ratio）
The scripts used in this paper to calculate dN/dS ratios are largely based on a previous analysis developed by Kistler and colleagues [1].

# Causal analysis
We conducted a nonlinear causality test based on CCM analysis using the 'rEDM' R package (Version 0.7.5) [2,3], which can be installed with the following R code: 

devtools::install_github("ha0ye/rEDM")


References
1. Kistler, K. E., Huddleston, J. & Bedford, T. Rapid and parallel adaptive mutations in spike S1 drive clade success in SARS-CoV-2. Cell Host Microbe 30, 545-555 e544 (2022).
2. Sugihara, G. et al. Detecting Causality in Complex Ecosystems. Science 338, 496-500 (2012).
3. Ye, H., Clark, A., Deyle, E., & Sugihara, G., rEDM: An R package for Empirical Dynamic Modeling and Convergent Cross Mapping. https://ha0ye.github.io/rEDM/articles/rEDM.html#generalized-takenss-theorem (2019). 
