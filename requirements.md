# Workflow requirements

## Dependencies installation

To install all the requisite packages used herein, please run the following:

```{r}
install.packages('BiocManager')
BiocManager::install(c("devtools", "lmtest", "parallel", "MAST", "MultiAssayExperiment",
                       "SummarizedExperiment", "tidyverse", "data.table",
                       "dplyr", "gplots", "ggplot2", "RColorBrewer",
                       "ggpubr", "tools", "splatter", "Matrix", "hypergate", "rio",
                       "car", "abind", "biomaRt", "scater", "gridExtra", "reshape",
                       "grid", "patchwork", "RSvgDevice", "scde", "UpSetR", "pROC",
                       "DESeq2", "tximport", "reshape2", "testthat", "truncnorm",
                       "MASS", "aggregation"))

## Overwrite install; Seurat must be <3.0.0
devtools::install_github('satijalab/seurat', ref = '65b77a9') # tagged release 2.3.3

devtools::install_github("thomasp85/patchwork") 

```

## Dependencies loading

To ensure all packages load, please run:

```{r}
library(lmtest)
library(parallel)
library(MAST)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)
library(dplyr)
library(patchwork)
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(tools)
library(splatter)
library(Matrix)
library(hypergate)
library(rio)
library(car)
library(abind)
library(biomaRt)
library(scater)
library(gridExtra)
library(reshape)
library(grid)
library(patchwork)
library(RSvgDevice)
library(Seurat)
library(scde)
library(UpSetR)
library(pROC)
library(DESeq2)
library(tximport)
library(reshape2)
library(testthat)
library(truncnorm)
library(MASS)
library(aggregation)
```

## Session Info

Code was tested on the following configuration:

```
R version 3.6.0 (2019-04-26)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 14.04.5 LTS

Matrix products: default
BLAS/LAPACK: /app/easybuild/software/OpenBLAS/0.2.18-GCC-5.4.0-2.26-LAPACK-3.6.1/lib/libopenblas_prescottp-r0.2.18.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      tools     stats4    parallel  stats     graphics  grDevices utils     datasets 
[10] methods   base     

other attached packages:
 [1] aggregation_1.0.1           MASS_7.3-51.4               truncnorm_1.0-8            
 [4] testthat_2.1.1              reshape2_1.4.3              tximport_1.12.0            
 [7] DESeq2_1.24.0               pROC_1.14.0                 UpSetR_1.3.3               
[10] scde_2.12.0                 flexmix_2.3-15              lattice_0.20-38            
[13] Seurat_2.3.3                cowplot_0.9.4               RSvgDevice_0.6.4.4         
[16] reshape_0.8.8               gridExtra_2.3               scater_1.12.0              
[19] biomaRt_2.40.0              abind_1.4-5                 car_3.0-2                  
[22] carData_3.0-2               rio_0.5.16                  hypergate_0.8.1            
[25] Matrix_1.2-17               splatter_1.8.0              ggpubr_0.2                 
[28] magrittr_1.5                RColorBrewer_1.1-2          gplots_3.0.1.1             
[31] patchwork_0.0.1             data.table_1.12.2           forcats_0.4.0              
[34] stringr_1.4.0               dplyr_0.8.1                 purrr_0.3.2                
[37] readr_1.3.1                 tidyr_0.8.3                 tibble_2.1.1               
[40] ggplot2_3.1.1               tidyverse_1.2.1             MultiAssayExperiment_1.10.0
[43] MAST_1.10.0                 SingleCellExperiment_1.6.0  SummarizedExperiment_1.14.0
[46] DelayedArray_0.10.0         BiocParallel_1.18.0         matrixStats_0.54.0         
[49] Biobase_2.44.0              GenomicRanges_1.36.0        GenomeInfoDb_1.20.0        
[52] IRanges_2.18.0              S4Vectors_0.22.0            BiocGenerics_0.30.0        
[55] lmtest_0.9-37               zoo_1.8-5                   here_0.1                   
[58] fs_1.3.1                    usethis_1.5.0               devtools_2.0.2             

loaded via a namespace (and not attached):
  [1] SparseM_1.77              prabclus_2.2-7            R.methodsS3_1.7.1        
  [4] acepack_1.4.1             bit64_0.9-7               knitr_1.22               
  [7] irlba_2.3.3               R.utils_2.8.0             Rook_1.1-1               
 [10] rpart_4.1-15              RCurl_1.95-4.12           generics_0.0.2           
 [13] metap_1.1                 snow_0.4-3                callr_3.2.0              
 [16] RSQLite_2.1.1             RANN_2.6.1                proxy_0.4-23             
 [19] bit_1.1-14                xml2_1.2.0                lubridate_1.7.4          
 [22] assertthat_0.2.1          viridis_0.5.1             xfun_0.7                 
 [25] RMTstat_0.3               hms_0.4.2                 DEoptimR_1.0-8           
 [28] progress_1.2.1            caTools_1.17.1.2          readxl_1.3.1             
 [31] geneplotter_1.62.0        igraph_1.2.4.1            DBI_1.0.0                
 [34] htmlwidgets_1.3           RcppArmadillo_0.9.400.3.0 backports_1.1.4          
 [37] trimcluster_0.1-2.1       annotate_1.62.0           gbRd_0.4-11              
 [40] quantreg_5.38             Cairo_1.5-10              remotes_2.0.4            
 [43] ROCR_1.0-7                withr_2.1.2               robustbase_0.93-5        
 [46] checkmate_1.9.3           prettyunits_1.0.2         mclust_5.4.3             
 [49] cluster_2.0.9             ape_5.3                   diffusionMap_1.1-0.1     
 [52] segmented_0.5-4.0         lazyeval_0.2.2            crayon_1.3.4             
 [55] genefilter_1.66.0         hdf5r_1.2.0               edgeR_3.26.1             
 [58] pkgconfig_2.0.2           nlme_3.1-140              vipor_0.4.5              
 [61] pkgload_1.0.2             blme_1.0-4                nnet_7.3-12              
 [64] rlang_0.3.4               diptest_0.75-7            MatrixModels_0.4-1       
 [67] extRemes_2.0-10           doSNOW_1.0.16             modelr_0.1.4             
 [70] rsvd_1.0.0                cellranger_1.1.0          rprojroot_1.3-2          
 [73] distillery_1.0-6          boot_1.3-22               base64enc_0.1-3          
 [76] beeswarm_0.2.3            ggridges_0.5.1            processx_3.3.1           
 [79] rjson_0.2.20              png_0.1-7                 viridisLite_0.3.0        
 [82] bitops_1.0-6              R.oo_1.22.0               Lmoments_1.3-1           
 [85] KernSmooth_2.23-15        blob_1.1.1                DelayedMatrixStats_1.6.0 
 [88] lars_1.2                  brew_1.0-6                scales_1.0.0             
 [91] memoise_1.1.0             plyr_1.8.4                ica_1.0-2                
 [94] bibtex_0.4.2              gdata_2.18.0              zlibbioc_1.30.0          
 [97] compiler_3.6.0            lsei_1.2-0                pcaMethods_1.76.0        
[100] lme4_1.1-21               fitdistrplus_1.0-14       cli_1.1.0                
[103] dtw_1.20-1                XVector_0.24.0            pbapply_1.4-0            
[106] ps_1.3.0                  htmlTable_1.13.1          Formula_1.2-3            
[109] mgcv_1.8-28               tidyselect_0.2.5          stringi_1.4.3            
[112] BiocSingular_1.0.0        locfit_1.5-9.1            latticeExtra_0.6-28      
[115] rstudioapi_0.10           foreach_1.4.4             foreign_0.8-71           
[118] scatterplot3d_0.3-41      Rtsne_0.15                digest_0.6.18            
[121] fpc_2.1-11.2              Rcpp_1.0.1                broom_0.5.2              
[124] SDMTools_1.1-221.1        httr_1.4.0                AnnotationDbi_1.46.0     
[127] npsurv_0.4-0              kernlab_0.9-27            Rdpack_0.11-0            
[130] colorspace_1.4-1          reticulate_1.12           rvest_0.3.3              
[133] XML_3.98-1.19             splines_3.6.0             sessioninfo_1.1.1        
[136] xtable_1.8-4              jsonlite_1.6              nloptr_1.2.1             
[139] modeltools_0.2-22         R6_2.4.0                  Hmisc_4.2-0              
[142] pillar_1.4.0              htmltools_0.3.6           glue_1.3.1               
[145] minqa_1.2.4               BiocNeighbors_1.2.0       class_7.3-15             
[148] codetools_0.2-16          tsne_0.1-3                pkgbuild_1.0.3           
[151] mvtnorm_1.0-10            mixtools_1.1.0            curl_3.3                 
[154] ggbeeswarm_0.6.0          gtools_3.8.1              zip_2.0.2                
[157] openxlsx_4.1.0            limma_3.40.0              survival_2.44-1.1        
[160] desc_1.2.0                munsell_0.5.0             GenomeInfoDbData_1.2.1   
[163] iterators_1.0.10          haven_2.1.0               gtable_0.3.0             
```
