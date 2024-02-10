# Global disparities in COVID-19 vaccine coverage associated with trajectories of SARS-CoV-2 adaptation

**Phylogenetic analysis**

We employed the NextStrain [1] pipeline (Version 6.2.0) to produce phylogenetic reconstructions for a specific country or region. **The phylogeny was constructed for the specific country using subsampled sequences and a small set of reference sequences, totaling 274 sequences maintained by the Nextstrain team as the genetic context. The sequences and metadata can be found at adaptation_calculation/ncov_Denmark/data/sequences.fasta and adaptation_calculation/ncov_Denmark/data/metadata.tsv.** We downloaded these sequences and metadata on April 2, 2023, from https://data.nextstrain.org/files/ncov/open/reference/metadata.tsv.xz and https://data.nextstrain.org/files/ncov/open/reference/sequences.fasta.xz.

Augur [2] (Version 21.1.0), IQ-TREE [3] (Version 2.1.2), and TreeTime [4] (Version 0.9.5) were used to infer time-resolved phylogenies. The joint inference mode in TreeTime was employed to infer ancestral sequences. We performed 10 phylogenetic reconstructions for each country based on the 10 bootstrap subsampling genome datasets. **As an example, we provided the phylogenetic tree for Denmark, which was constructed using one of the 10 bootstrap subsampling genome datasets for Denmark. Note that if you want to visualize the phylogenetic tree, unzip adaptation_calculation/ncov_Denmark/auspice/ncov_Denmark-build.zip and drop its contents onto https://auspice.us/.**

Given that the strategy of phylogeny reconstruction may impact the robustness of results, we collected a set of publicly available SARS-CoV-2 phylogenetic analyses for specific countries on January 2, 2024. These SARS-CoV-2 phylogeny are produced and maintained by the Nextstrain groups and include countries from Africa (https://nextstrain.org/groups/africa-cdc, created as part of [5]), Europe and North America (https://nextstrain.org/groups/neherlab). In total, 50 countries align with the scope of our investigation, with the source of website links available in Table S1 of the supplemental information. We provide the phylogenetic tree for Denmark, downloaded on 2 January 2024, as an example (adaptation_calculation/ncov_neherlab_Denmark/auspice/neherlab_ncov_denmark.json). If you want to visualize it, drop it onto https://auspice.us/

**Calculation of the ratio of nonsynonymous to synonymous divergence （dN/dS ratio）,mutation accumulation and mutational fitness**

We calculated the dN/dS ratios to assess selection pressure on SARS-CoV-2 proteins per country over time. For a specific country, the ratio was calculated by the accumulating nonsynonymous divergence sites relative to synonymous divergence sites between internal branches inferred from the specific country and the tree root in the phylogeny, using a time window of 2 months with a 1-month overlap. Furthermore, we measured the number of nonsynonymous mutations (by tallying the accumulation between the phylogeny root and internal branches) and the metric of mutational fitness to suggest SARS-CoV-2 adaptation. Mutations were categorized for each gene according to the Wuhan-Hu-1 reference sequence (found in “adaptation_calculation/ncov_Denmark/analysis/reference_seq_edited.gb”) to identify whether they were synonymous or nonsynonymous. The scripts used in this study for calculating dN/dS ratios and quantifying the number of nonsynonymous mutations largely rely on previous methods developed by Kistler and colleagues [6]. Mutational fitness was performed using Nextstrain pipelines (Version 6.2.0), estimated based on hierarchical Bayesian multinomial logistic regression that calculates the relative impact of individual mutations on the growth advantage of SARS-CoV-2 lineages [7].

**We provided the code for calculating adaptation metrics based on the phylogenetic analysis of Denmark to calculate dN/dS ratio (adaptation_calculation/ncov_Denmark/analysis/dnds_nextstrain_Denmark.ipynb), the number of nonsynonymous mutations, and mutational fitness (adaptation_calculation/ncov_Denmark/analysis/mutation_nextstrain_Denmark.ipynb).** The calculation of the three metrics of SARS-CoV-2 adaptation for other countries and additional subsampled datasets from the 10 bootstrap is consistent with the code used for Denmark, with the only difference being the read of the phylogenetic data.

**Additionally, we provided the code for calculating the dN/dS ratio (located in adaptation_calculation/ncov_neherlab_Denmark/analysis/dnds_nextstrain_Denmark.ipynb), the number of nonsynonymous mutations and mutational fitness (found in adaptation_calculation/ncov_neherlab_Denmark/analysis/mutation_nextstrain_Denmark.ipynb) based on the publicly available SARS-CoV-2 phylogeny.**

**Statistical analysis**

The distributed lag non-linear model (DLNM) from the “dlnm” R package (Version 2.4.7) [8] was fitted to investigate the potential non-linear association between selection pressure on SARS-CoV-2 regions and immunity from COVID-19 vaccines and natural infections.

**Causal analysis**

We conducted a nonlinear causality test based on convergent cross-mapping (CCM) analysis using the 'rEDM' R package (Version 0.7.5) [9,10], which can be installed with the following R code: 

devtools::install_github("ha0ye/rEDM")


**References**
1.	Hadfield, J. et al. Nextstrain: real-time tracking of pathogen evolution. Bioinform 34, 4121-4123 (2018).
2.	Huddleston, J. et al. Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens. J. Open Source Softw 6, 2906 (2021).
3.	Nguyen, L. T., Schmidt, H. A., von Haeseler, A. & Minh, B. Q. IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Mol Biol Evol 32, 268-274 (2015).
4.	Sagulenko, P., Puller, V. & Neher, R. A. TreeTime: Maximum-likelihood phylodynamic analysis. Virus Evol 4, vex042 (2018).
5.	Tegally, H. et al. The evolving SARS-CoV-2 epidemic in Africa: Insights from rapidly expanding genomic surveillance. Science 378, eabq5358 (2022).
6.  Kistler, K. E., Huddleston, J. & Bedford, T. Rapid and parallel adaptive mutations in spike S1 drive clade success in SARS-CoV-2. Cell Host Microbe 30, 545-555 e544 (2022).
7.	Obermeyer, F. et al. Analysis of 6.4 million SARS-CoV-2 genomes identifies mutations associated with fitness. Science 376, 1327-1332 (2022).
8.	Gasparrini, A. Modeling exposure-lag-response associations with distributed lag non-linear models. Stat Med 33, 881-899 (2014).
9. Sugihara, G. et al. Detecting Causality in Complex Ecosystems. Science 338, 496-500 (2012).
10. Ye, H., Clark, A., Deyle, E., & Sugihara, G., rEDM: An R package for Empirical Dynamic Modeling and Convergent Cross Mapping. https://ha0ye.github.io/rEDM/articles/rEDM.html#generalized-takenss-theorem (2019). 
