## Use Seurat 2.3.4
## source("http://bit.ly/archived-seurat")

#Load and parse data ----------------------------------------------------
library(MultiAssayExperiment)
## data <- readRDS('/home/lynnyi/NYMP_2018/embryo/EMTAB3929.rds')
data <- readRDS('./EMTAB3929.rds') ##EB: Downloaded from http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/EMTAB3929.rds
gene <- experiments(data)[['gene']]
tpms <- assays(gene)[['TPM']][,1:271]
counts <- assays(gene)[['count']][,1:271]
day <- c(rep(3, 81), rep(4, 190))
cells.1=which(day==3)
cells.2=which(day==4)

## Fixing the issue with subsetted genes
#Load and parse data ----------------------------------------------------
library(MultiAssayExperiment)
## data <- readRDS('/home/lynnyi/NYMP_2018/embryo/EMTAB3929.rds')

## Needs Seurat 2.3.4 not 3.0
# Use Seurat to call DESeq2 ---------------------------------------------
library(Seurat)
nbt = CreateSeuratObject(raw.data=round(counts), min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(counts) 
## deseq2 <- FindMarkers(nbt, ident.1=3, test.use='DESeq2', logfc.threshold=-Inf, min.pct=-Inf)
## saveRDS(deseq2, './deseq2_fixed.rds')
deseq2=DESeq2DETest(nbt,cells.1=cells.1,cells.2=cells.2)
saveRDS(deseq2, './deseq2_fixed2.rds')

# Use Seurat to call MAST -----------------------------------------------
nbt = CreateSeuratObject(raw.data=log(tpms+1), min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(tpms)
mast <- FindMarkers(nbt, ident.1=3, test.use='MAST', logfc.threshold=-Inf, min.pct=-Inf)
saveRDS(mast, './mast_fixed.rds')

# Use Seurat to run Tobit model --------------------------------------------
nbt = CreateSeuratObject(raw.data=tpms, min.cells=0, min.genes=0, project='Embryo')
nbt@ident <- as.factor(day)
names(nbt@ident) <- colnames(tpms)
##tobit <- FindMarkers(nbt, ident.1=3, test.use='tobit', logfc.threshold=-Inf, min.pct=-Inf)
tobit=TobitTest(nbt, cells.1 = cells.1, cells.2 = cells.2)
saveRDS(tobit, './tobit_fixed.rds')
