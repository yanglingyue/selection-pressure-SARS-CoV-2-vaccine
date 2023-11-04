# Global disparities in COVID-19 vaccine coverage associated with trajectories of SARS-CoV-2 adaptation

**Phylogenetic analysis**

We employed the NextStrain [1] pipeline (Version 6.2.0) to produce phylogenetic reconstructions for a specific country or region. The phylogeny was constructed for a specific country with subsampled sequences and some sequences around the world to place the country’s sequences under the context of the global pandemic. We downloaded sequences maintained by the Nextstrain team as the genetic context, which can be found at https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html. Augur [2] (Version 21.1.0), IQ-TREE [3] (Version 2.1.2), and TreeTime [4] (Version 0.9.5) were used to infer time-resolved phylogenies. The joint inference mode in TreeTime was employed to infer ancestral sequences.

**Calculation of the ratio of nonsynonymous to synonymous divergence （dN/dS ratio）**

The scripts used in this paper to calculate dN/dS ratios are largely based on a previous analysis developed by Kistler and colleagues [5].

**Statistical analysis**

The distributed lag non-linear model (DLNM) from the “dlnm” R package (Version 2.4.7) [6] was fitted to investigate the potential non-linear association between selection pressure on SARS-CoV-2 regions and immunity from COVID-19 vaccines and natural infections.

**Causal analysis**

We conducted a nonlinear causality test based on convergent cross-mapping (CCM) analysis using the 'rEDM' R package (Version 0.7.5) [6,7], which can be installed with the following R code: 

devtools::install_github("ha0ye/rEDM")


**References**
1.	Hadfield, J. et al. Nextstrain: real-time tracking of pathogen evolution. Bioinform 34, 4121-4123 (2018).
2.	Huddleston, J. et al. Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. J. Open Source Softw 6, 2906 (2021).
3.	Nguyen, L. T., Schmidt, H. A., von Haeseler, A. & Minh, B. Q. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 32, 268-274 (2015).
4.	Sagulenko, P., Puller, V. & Neher, R. A. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evol 4, vex042 (2018).
5.  Kistler, K. E., Huddleston, J. & Bedford, T. Rapid and parallel adaptive mutations in spike S1 drive clade success in SARS-CoV-2. Cell Host Microbe 30, 545-555 e544 (2022).
6.	Gasparrini, A. Modeling exposure-lag-response associations with distributed lag non-linear models. Stat Med 33, 881-899 (2014).
7. Sugihara, G. et al. Detecting Causality in Complex Ecosystems. Science 338, 496-500 (2012).
8. Ye, H., Clark, A., Deyle, E., & Sugihara, G., rEDM: An R package for Empirical Dynamic Modeling and Convergent Cross Mapping. https://ha0ye.github.io/rEDM/articles/rEDM.html#generalized-takenss-theorem (2019). 
