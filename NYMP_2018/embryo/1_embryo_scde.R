# Perform SCDE on embryo dataset -----------------------------
library(scde)
setwd("/fh/fast/gottardo_r/ebecht_working/logistic_regression_public/NYMP_2018/embryo/")

options(mc.cores = 1L)

genetpms <- readRDS("./embryonic_gene_tpms.rds")
coldata <- readRDS("./embryonic_metadata.rds")
## sg <- factor(coldata$day) names(sg) <- as.character(coldata$sample)
sg <- setNames(factor(coldata), colnames(genetpms))

genetpms <- apply(genetpms, 2, function(x) {
    storage.mode(x) <- "integer"
    x
})
colnames(genetpms) <- names(sg)

write_dir <- "."

## Requires flexmix 2.3-13 (above doesn't work)
error_models <- scde.error.models(counts = genetpms, groups = sg, n.cores = 4, threshold.segmentation = TRUE, 
    save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
saveRDS(error_models, file.path(write_dir, "scde_error_models.rds"))

prior <- scde.expression.prior(models = error_models, genetpms, show.plot = FALSE)
scde_DE <- scde.expression.difference(error_models, genetpms, prior, sg, verbose = 1, 
    n.cores = 4, n.randomizations = 100)
scde_DE$pvalue <- 2 * pnorm(-abs(scde_DE$Z))
saveRDS(scde_DE, file.path(write_dir, "scde.rds"))
