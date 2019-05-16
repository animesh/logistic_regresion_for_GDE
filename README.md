This repository contains code to reproduce our comparison of MAST-Sidak and multivariate logistic regression (reply to Ntranos, Yi et al, Nature Methods, 2019).

Every script should be run with R's current working directory inside the directory where the script is located. This is performed ad-hoc by the `setwd()` command, assuming that R begins at the top of the directory.

R package requirements, including a convenience wrapper for package install and loading, as well as session info for the environment in which the code was tested, are enumerated in the `requirements.md` file.


The order in which to recapitulate the manuscript's analysis is as follows:

1. Run the reproduced and extended figure about CD45:
./NYMP_2018/10x_example-logR/script_EB.R

2. Run data processing scripts:
./figshare/EB/simulations_analysis_clean.R (simulated datasets)
./figshare/EB/shuffling_analysis.R (empirical datasets)


3. Run the adapted analysis from Ntranos et al for the analysis of the embryo dataset
./NYMP_2018/embryo/0_embryo_seurat.R
./NYMP_2018/embryo/1_embryo_scde.R
./NYMP_2018/embryo/2_embryo_LR.R
./NYMP_2018/embryo/3_plot_figures_EB.R 

4. Finally, run the script corresponding to figures from the paper:
./paper_v1/scripts_v2.R

