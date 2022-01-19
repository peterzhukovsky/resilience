This repository contains a set of scripts used to analyze resilience to amyliod-beta accumulation and resilience conferred by education in OASIS-3 cross-sectional data.

# 1. Functional connectivity signatures of resilience conferred by education
*analysis_ICA25d.m* includes scripts that import the data, exclude several participants identified by QC; and run linear models (*fitlm*) to test for the effects of education on pairwise functional connectivity. The functional connectivity variables undergo exclusion of outliers (>3SD from the mean) first to ensure data quality is adequate. 

Next , we use 1,000 permutations to test for significance of the effect of education. Each permutation involves reshuffling the outcome variable (of the 210 functional connectivities). From each permuted model we get a "by chance" F-value which allows us to generate a permutation distribution. 

We provide code for visualizing the connectivity matrices and circular graphs. 

# 2. Functional connectivities predicting memory recall
*pls_ICA25d.m* set of scripts includes code that imports and cleans up data similar to *ICA25d.m*. USing the 210 connectivities as predictors (X), and age and delayed memory recall as outcomes (Y), partial least squares PLS regression is ran. Permutation testing n=5,000 is used to test whether the PLS explains a significant amount of variance in the outcome (Y) and bootstrapping is used to test which connectivities show a robust contribution (Z score > 3) to the latent PLS variable.

# 3. Functional connectivities mediate the relationship between mean cortical amyloid-beta and memory recall
A PLS structural equation model script is used here (*pls_sem_O3.m*). PLS-SEM toolbox outputs are provided in the folder called *report_FINAL_BL* for the baseline data and *report_extended_sample* for the extension analysis.
https://www.mathworks.com/matlabcentral/fileexchange/54147-pls-sem-toolbox

# 4. Extension analyses
Scripts used to run the analyses in points 2. and 3. (i.e. the PLS and PLS-SEM) in the extension sample are called *pls_ICA25d_m.m* and *pls_sem_O3.m*, respectively. 
