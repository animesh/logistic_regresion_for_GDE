## EB
t2g <- readRDS("./t2g.rds")
## 

# Load methods' results ----------------------------------------------
scde <- readRDS("./scde.rds")
mast <- readRDS("./mast.rds")
tobit <- readRDS("./tobit.rds")
deseq2 <- readRDS("./deseq2.rds")
LR <- readRDS("./LR.rds")


## EB
LR$genes <- as.character(LR$genes)
egv2eg <- function(genes) {
    sapply(strsplit(genes, ".", fixed = TRUE), "[", 1)
}
rownames(scde) <- egv2eg(rownames(scde))
rownames(mast) <- egv2eg(rownames(mast))
rownames(tobit) <- egv2eg(rownames(tobit))
rownames(deseq2) <- egv2eg(rownames(deseq2))
rownames(LR) <- LR$genes
## /EB

## EB
names_pretty <- setNames(c("log reg", "SCDE", "monocle - Tobit", "DESeq2", "MAST"), 
    c("LR", "scde", "tobit", "deseq2", "mast"))
before_subsetting <- sapply(list("scde", "mast", "tobit", "deseq2", "LR"), function(x) {
    setNames(nrow(get(x, envir = .GlobalEnv)), names_pretty[x])
})
## /EB

# Make sure the genes match, since some pulled from conquer and some genes (LR)
# pulled from ensembl
scde <- scde[rownames(scde) %in% LR$genes, ]
tobit <- tobit[rownames(tobit) %in% LR$genes, ]
deseq2 <- deseq2[rownames(deseq2) %in% LR$genes, ]
mast <- mast[rownames(mast) %in% LR$genes, ]

# Perform multiple hypothesis adjustment ----------------------------
LR$p_val_adj <- p.adjust(LR$pvalues, "BH")
scde$p_val_adj <- p.adjust(scde$pvalue, "BH")

# Create upset plot Supp Fig 6 -------------------------------- Take top 3000
# genes from all methods
l <- list(LR, scde, tobit, deseq2, mast)
## EB
names(l) <- c("log reg", "SCDE", "monocle - Tobit", "DESeq2", "MAST")
## /EB
test <- lapply(l, function(x) rownames(x[order(x$p_val_adj), ])[1:3000])
names(test) <- c("log reg", "SCDE", "monocle - Tobit", "DESeq2", "MAST")
library(UpSetR)
tiff("./upset.tiff", height = 3000, width = 4000, res = 500)
upset(fromList(test), order.by = "freq", text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()

## EB : Number of features assayed in the plot for each method

# Find what is unique to logistic regression
unique <- setdiff(test[[1]], test[[2]])
unique <- setdiff(unique, test[[3]])
unique <- setdiff(unique, test[[4]])
unique <- setdiff(unique, test[[5]])

library(reshape)
library(ggplot2)
plot_gene <- function(gene_name) {
    i <- which(t2g$ensembl_gene_id == gene_name)
    tx <- t2g$ensembl_transcript_id[i]
    print(length(tx))
    tx_counts <- tx_counts[rownames(tx_counts) %in% tx, ]
    tx_counts <- log(tx_counts + 1)
    tx_counts <- melt(tx_counts)
    names(tx_counts) <- c("transcript", "day", "counts")
    tx_counts$day <- as.factor(ifelse(grepl("E3", tx_counts$day), "E3", "E4"))
    g <- ggplot(tx_counts, aes(transcript, counts, fill = day)) + geom_boxplot(alpha = 0.2) + 
        ylab("log(counts+1)")
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g <- g + ggtitle(gene_name)
    tiff(paste0(gene_name, ".tiff"))
    print(g)
    dev.off()
}

plot_gene_violin <- function(gene_name) {
    i <- which(t2g$ensembl_gene_id == gene_name)
    tx <- t2g$ensembl_transcript_id[i]
    print(length(tx))
    tx_counts <- tx_counts[rownames(tx_counts) %in% tx, ]
    tx_counts <- log(tx_counts + 1)
    tx_counts <- melt(tx_counts)
    names(tx_counts) <- c("transcript", "day", "counts")
    tx_counts$day <- as.factor(ifelse(grepl("E3", tx_counts$day), "E3", "E4"))
    g <- ggplot(tx_counts, aes(transcript, counts, fill = day)) + geom_violin(alpha = 0.2) + 
        ylab("log(counts+1)")
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    g <- g + ggtitle(gene_name)
    tiff(paste0("violin_", gene_name, ".tiff"))
    print(g)
    dev.off()
}

# Supp Fig 6b
plot_gene("ENSG00000204481")
plot_gene_violin("ENSG00000204481")

# Supp Fig 6c
plot_gene("ENSG00000128253")
plot_gene_violin("ENSG00000128253")

# Supp Fig 6d
plot_gene("ENSG00000122565")
plot_gene_violin("ENSG00000122565")

# Supp Fig 6e
plot_gene("ENSG00000106682")
plot_gene_violin("ENSG00000106682")

## EB
shared <- rownames(l[[1]])
for (i in 2:length(l)) {
    shared <- intersect(shared, rownames(l[[i]]))
}
test2 <- lapply(l, function(x) {
    x <- x[shared, ]
    rownames(x[order(x$p_val_adj), ])[1:3000]
})
names(test2) <- c("log reg", "SCDE", "monocle - Tobit", "DESeq2", "MAST")
library(UpSetR)
tiff("./upset_EB.tiff", height = 3000, width = 4000, res = 500)
upset(fromList(test2), order.by = "freq", text.scale = c(2, 2, 2, 1.3, 2, 1.2))
dev.off()

before_subsetting <- sapply(list("scde", "mast", "tobit", "deseq2", "LR"), function(x) {
    setNames(nrow(get(x, envir = .GlobalEnv)), names_pretty[x])
})
after_subsetting <- sapply(l, function(x) {
    nrow(x[shared, ])
})

png("../../paper_v1/graphs/Sup fig 1A.png", width = 2000, res = 300, height = 2000)
par(bty = "l", mar = c(8, 4, 2, 2))
barplot(before_subsetting, ylab = "Number of features assayed", las = 3)
dev.off()

png("../../paper_v1/graphs/Sup fig 1B.png", width = 2000, res = 300, height = 2000)
par(bty = "l", mar = c(8, 4, 2, 2))
barplot(after_subsetting, ylab = "Number of features assayed", las = 3)
dev.off()

## EB Number of 0s per tx per condition
w <- grepl("E4", colnames(tx_counts))
E4_0 <- (tx_counts[, w] == 0)
E3_0 <- (tx_counts[, !w] == 0)

E4_ag <- aggregate_df(E4_0, t2g[match(rownames(E4_0), t2g[, 2]), 1], sum)
E3_ag <- aggregate_df(E3_0, t2g[match(rownames(E3_0), t2g[, 2]), 1], sum)

E4_ag <- rowSums(E4_ag)
E3_ag <- rowSums(E3_ag)

wmast <- rownames(subset(l[[5]][shared, ], p_val < 0.05))
wlr <- rownames(subset(l[[1]][shared, ], pvalues < 0.05))

png(width = 2000, height = 2000, res = 300, "counting_0s.png")
par(mfrow = c(2, 2))
freqplot(E4_ag[intersect(wmast, wlr)], E3_ag[intersect(wmast, wlr)])
title(main = "D.E. for both LR & MAST")
abline(a = 0, b = sum(!w)/sum(w))
title(xlab = "Sum of 0 in tx counts in E4 cells [ag to gene]", ylab = "Sum of 0 in tx counts in E3 cells [ag to gene]")
freqplot(E4_ag[setdiff(wmast, wlr)], E3_ag[setdiff(wmast, wlr)])
title(main = "Unique to MAST")
abline(a = 0, b = sum(!w)/sum(w))
freqplot(E4_ag[setdiff(wlr, wmast)], E3_ag[setdiff(wlr, wmast)])
title(main = "Unique to LR")
abline(a = 0, b = sum(!w)/sum(w))
freqplot(E4_ag[setdiff(shared, c(wmast, wlr))], E3_ag[setdiff(shared, c(wmast, wlr))])
title(main = "Negative for both")
abline(a = 0, b = sum(!w)/sum(w))
dev.off()
