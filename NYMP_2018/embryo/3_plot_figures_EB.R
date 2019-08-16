####################
## Reproduction and modification of Supplementary Figure 6 of the Nnatros, Yi, Melsted and Pachter 2019 paper.
## This script is adapted from https://github.com/pachterlab/NYMP_2018/blob/master/embryo/3_plot_figures.R
## It's most important modification is to appropriately compute intersections on a consistent set of input features
## Modifications from the original script are listed using ## EB and ## /EB flags
## You should run the rest (scripts 0, 1 and 2) of the https://github.com/pachterlab/NYMP_2018/tree/master/embryo directory before running this script
####################

## Transcripts to gene mapping
t2g=readRDS("./t2g.rds")

# Load methods' results ----------------------------------------------
scde <- readRDS('./scde.rds')
mast <- readRDS('./mast.rds')
tobit <- readRDS('./tobit.rds')
deseq2 <- readRDS('./deseq2.rds')
LR <- readRDS('./LR.rds')


##EB: Making sure we use the same format for gene IDs (some are ensembl_ids.version and some are only ensembl_ids)
LR$genes=as.character(LR$genes)
egv2eg=function(genes){
    sapply(strsplit(genes,".",fixed=TRUE),"[",1)
}
rownames(scde)=egv2eg(rownames(scde))
rownames(mast)=egv2eg(rownames(mast))
rownames(tobit)=egv2eg(rownames(tobit))
rownames(deseq2)=egv2eg(rownames(deseq2))
rownames(LR)=LR$genes
## /EB

## EB: Rewriting methods' names
names_pretty=setNames(c('logistic\nregression', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST'),c("LR","scde","tobit","deseq2","mast"))
before_subsetting=sapply(list("scde","mast","tobit","deseq2","LR"),function(x){setNames(nrow(get(x,envir=.GlobalEnv)),names_pretty[x])})
## /EB

# Make sure the genes match, since some pulled from conquer and some genes (LR) pulled from ensembl 
scde <- scde[rownames(scde) %in% LR$genes,]
tobit <- tobit[rownames(tobit) %in% LR$genes,]
deseq2 <- deseq2[rownames(deseq2) %in% LR$genes,]
mast <- mast[rownames(mast) %in% LR$genes,]

# Perform multiple hypothesis adjustment ----------------------------
LR$p_val_adj <- p.adjust(LR$pvalues, 'BH')
scde$p_val_adj <- p.adjust(scde$pvalue, 'BH')

##EB
library(RSvgDevice)
devSVG("../../paper_v2/graphs/1C_pretty.svg",width=3.5,height=5)
par(bty="l",mar=c(16,8,2,2),cex.lab=1.5)
barplot(before_subsetting[c("logistic\nregression","MAST","SCDE","monocle - Tobit","DESeq2")],ylab="Number of input genes",las=3,space=0.1,cex.names=1.5,cex.axis=1.2,width=.3)
dev.off()
##/EB

# Create upset plot Supp Fig 6 --------------------------------
# Take top 3000 genes from all methods
l <- list(LR, scde, tobit, deseq2, mast)
## EB
names(l)=c('log reg', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
## /EB
test <- lapply(l, function(x) rownames(x[order(x$p_val_adj),])[1:3000])
names(test) <- c('log reg', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
library(UpSetR)
tiff('./upset.tiff', height=3000, width = 4000, res= 500)
upset(fromList(test), order.by='freq', text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()
## EB : Number of features assayed in the plot for each method

paper_genes=c("ENSG00000204481","ENSG00000128253","ENSG00000122565","ENSG00000106682")

####################
## Updated figures after removing the  pre-filtering
####################
## EB
t2g=readRDS("./t2g.rds")
## /EB

# Load methods' results ----------------------------------------------
scde <- readRDS('./scde.rds')
mast <- readRDS('./mast_fixed.rds')
tobit <- readRDS('./tobit_fixed.rds')
deseq2 <- readRDS('./deseq2_fixed2.rds')
LR <- readRDS('./LR.rds')

##EB: Making sure we use the same format for gene IDs (some are ensembl_ids.version and some are only ensembl_ids)
LR$genes=as.character(LR$genes)
egv2eg=function(genes){
    sapply(strsplit(genes,".",fixed=TRUE),"[",1)
}
rownames(scde)=egv2eg(rownames(scde))
rownames(mast)=egv2eg(rownames(mast))
rownames(tobit)=egv2eg(rownames(tobit))
rownames(deseq2)=egv2eg(rownames(deseq2))
rownames(LR)=LR$genes
## /EB

## EB: Rewriting methods' names
names_pretty=setNames(c('logistic\nregression', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST'),c("LR","scde","tobit","deseq2","mast"))
before_subsetting=sapply(list("scde","mast","tobit","deseq2","LR"),function(x){setNames(nrow(get(x,envir=.GlobalEnv)),names_pretty[x])})
## /EB

# Perform multiple hypothesis adjustment ----------------------------
LR$p_val_adj <- p.adjust(LR$pvalues, 'BH')
scde$p_val_adj <- p.adjust(scde$pvalue, 'BH')
deseq2$p_val_adj <- p.adjust(deseq2$p_val, 'BH')
tobit$p_val_adj <- p.adjust(tobit$p_val, 'BH')

# Create upset plot Supp Fig 6 --------------------------------
# Take top 3000 genes from all methods
l <- list(LR, scde, tobit, deseq2, mast)
## EB
names(l)=c('log reg', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
## /EB

library(reshape)
library(ggplot2)

## EB
## Plotting updated Fig S6
shared=rownames(l[[1]])
for(i in 2:length(l)){
    shared=intersect(shared,rownames(l[[i]]))
}
library(UpSetR)
library(RSvgDevice)

before_subsetting=sapply(list("scde","mast","tobit","deseq2","LR"),function(x){setNames(nrow(get(x,envir=.GlobalEnv)),names_pretty[x])})
after_subsetting=sapply(l,function(x){nrow(x[shared,])})

test2 <- lapply(l, function(x){
    x=x[shared,]
    rownames(x[order(x$p_val_adj),])[1:3000]
})

names(test2) <- c('logistic\nregression', 'SCDE', 'monocle - Tobit', 'DESeq2', 'MAST')
devSVG('../../paper_v2/graphs/1D_UpSet.svg')
upset(fromList(test2), order.by='freq', text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()

sup_table_1=sapply(l,function(x){x[paper_genes,"p_val_adj"]})
rownames(sup_table_1)=paper_genes
