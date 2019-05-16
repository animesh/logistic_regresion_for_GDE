## =============================================================================
## This script reproduces and extends the histograms of Figure 2 of the Ntranos,
## Yi, Melsted and Pachter (Nature Methods, 2019) paper. It is adapted and translated
## from the python script available at:
## https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
##
## This script performs the following:
##   1. Data import
##   2. Train different logistic regressions (univariate on genes, univariate
##      on TCCs with Bonferonni aggregation to the gene leve, multivariate on
##      TCCs (aggregation to genes using likelihood ratio tests). In addition
##      to the original paper, we also apply MAST on TCCs followed by Sidak
##      aggregation of the p-values to the gene level. These models are trained
##      on 200 subsamples of sizes 2,000 (out of 3,000) cells per cell type,
##      for each pair of cell type.
##   3. Plot the distributions of p-values as histograms
##
## =============================================================================

library(Matrix)
library(tidyverse)
library(MultiAssayExperiment)
library(MAST)
library(data.table)
library(biomaRt)
library(tools)

source("wrappers.R")
setwd("NYMP_2018/10x_example-logR/")
options("mc.cores" = parallel::detectCores())


## =============================================================================
## 1. Data import and reshape
## =============================================================================

## -----------------------------------------------------------------------------
## Load data.
## Source: https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/matrix.tcc.mtx.gz
X <- readMM("./matrix.tcc.mtx")
colnames(X) <- 1:ncol(X) - 1

## Creating mappings
## Transcripts to ensembl transcript IDs
TX_TO_ENST <- read.csv("./TX_to_ENST.csv", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
TX_TO_ENST <- setNames(TX_TO_ENST[, 2], TX_TO_ENST[, 1])

## Transcripts to ensembl gene IDs
TX_TO_ENSG <- read.csv("./TX_to_ENSG.csv", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
TX_TO_ENSG <- setNames(TX_TO_ENSG[, 2], TX_TO_ENSG[, 1])

## Ensembl gene IDs to HUGO symbols
ENSG_to_name <- read.csv("./ENSG_to_name.csv", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
ENSG_to_name <- setNames(ENSG_to_name[, 2], ENSG_to_name[, 1])

## TCC to TX mappings
EC_dict <- read.csv("./matrix.ec", sep = "\t", stringsAsFactors = FALSE, header = FALSE)
EC_dict <- setNames(strsplit(EC_dict[, 2], ","), EC_dict[, 1])
EC_to_tx <- stack(EC_dict, stringsAsFactors = FALSE)
EC_to_tx <- as.matrix(EC_to_tx)
mode(EC_to_tx) <- "character"
colnames(EC_to_tx) <- c("TX", "TCC")

## TCC to HUGO symbol mapping
EC_to_gene_names <- stack(EC_dict, stringsAsFactors = FALSE)
colnames(EC_to_gene_names) <- c("TX", "TCC")
EC_to_gene_names <- apply(EC_to_gene_names, 2, as.character)
EC_to_gene_names <- cbind(EC_to_gene_names, ENSG = TX_TO_ENSG[EC_to_gene_names[, "TX"]])
EC_to_gene_names <- cbind(EC_to_gene_names, name = ENSG_to_name[EC_to_gene_names[, "ENSG"]])

## HUGO to TCCs
gene_names_to_ECs <- split(EC_to_gene_names[, "TCC"], EC_to_gene_names[, "name"])
gene_names_to_ECs <- lapply(gene_names_to_ECs, unique)
gene_names_to_ECs <- lapply(gene_names_to_ECs, sort)

## TCC to HUGO symbol mapping (finishing)
EC_to_gene_names <- split(EC_to_gene_names[, "name"], EC_to_gene_names[, "TCC"])
EC_to_gene_names <- lapply(EC_to_gene_names, unique)
EC_to_gene_names <- lapply(EC_to_gene_names, sort)

## -----------------------------------------------------------------------------
## Aggregating TCC matrix to transcripts matrix (Takes a while)
file <- "matrix_transcripts.rds" ## load precomputed results
if (file.exists(file)) {
  X_tx <- readRDS(file)
} else {
  X_mat <- as.matrix(X)
  X_split <- split_matrix(X_mat, colnames(X_mat), byrow = FALSE) ## split the matrix
  i <- 0
  ## For each TCC, create a matrix whose rows are the transcripts mapping to that
  ## TCC, and columns are single cell. Count is the count vector of TCC, divided
  ## by the number of transcripts, and repeated for each transcript (so we allocate
  ## 1/n count to each transcript where n is the number of transcripts this TCC maps to)
  X_split_tx <- lapply(
    X_split,
    function(x) {
      i <<- i + 1
      if (i %% 1000 == 0) {
        print(i)
      }
      ec <- colnames(x)
      txs <- EC_dict[[ec]]
      xp <- x[, 1] / length(txs)
      matrix(nrow = nrow(x), ncol = length(txs), data = rep(xp, length(txs)), byrow = FALSE, dimnames = list(rownames(x), txs))
    }
  )
  tx_split <- do.call(c, lapply(X_split_tx, colnames)) ## Tabulate for each TCC the mapping transcripts
  X_split_tx <- do.call(cbind, X_split_tx) ## Concatenate matrix, columns are now transcripts and rows are single cells. Transcripts who are mapped by distinct TCCs will appear multiple times
  X_split_tx <- split_matrix(X_split_tx, tx_split, byrow = FALSE) ## Split matrix to a list of matrices mapping to a given transcript
  X_split_tx <- lapply(X_split_tx, rowSums) ## Aggregate matrices for each transcript (sum counts)
  X_split_tx <- do.call(rbind, X_split_tx) ## Concatenate matrix
  rownames(X_split_tx) <- TX_TO_ENST[rownames(X_split_tx)] ## Set ensembl transcript IDs as rownames
  X_tx <- t(Matrix(X_split_tx))
  rm(X_split_tx)
  saveRDS(X_tx, file = file)
}

## Find transcripts and gene ensembl IDs associated with CD45 (PTPRC gene)
## Result is: "ENST00000348564", "ENST00000442510", "ENST00000530727", "ENST00000367367",
## "ENST00000413409", "ENST00000529828", "ENST00000575923", "ENST00000573679",
## "ENST00000576833", "ENST00000573477", "ENST00000571847", "ENST00000574441",
## as expected per:
## https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
mart <- useMart("ensembl")
mart <- useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
cd45_ensg <- names(ENSG_to_name[ENSG_to_name == "PTPRC"])
annot <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
               filters = "ensembl_gene_id", mart = mart, values = cd45_ensg)
annot <- setNames(annot[, "ensembl_gene_id"], annot[, "ensembl_transcript_id"])
cd45_txs <- intersect(names(annot), colnames(X_tx))
cd45_txs <- cd45_txs[Matrix:::colSums(X_tx[, cd45_txs]) > 0.003 * nrow(X_tx)]


## -----------------------------------------------------------------------------
## Filtering on expressed genes (>0.003 counts per cell) and uniquely-mapping TCCs (gene-wise)
## Making sure we have the right counts for each TCC mapping to PTPRC.
## See code block 12 in source:
## https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
cbind(gene_names_to_ECs[["PTPRC"]], Matrix:::colSums(X[, gene_names_to_ECs[["PTPRC"]]])) 

## Plotting function: filters multimapping TCCs for a gene, plots counts and
## threshold of mean counts per cell, returns vector of TCCs passing both
## the non-multimapping and minimum mean count per cell criteria
## See code block 13 of source below for the same result in Python by Ntranos, Yi et al.:
## https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb 
plot_gene_filt <- function(gene, threshold = 0.003) {
  ec_counts <- Matrix:::colSums(X[, gene_names_to_ECs[[gene]]])
  names(ec_counts) <- gene_names_to_ECs[[gene]]
  ec_multimapping <- sapply(EC_to_gene_names[as.character(gene_names_to_ECs[[gene]])], length) > 1

  o <- order(ec_counts, decreasing = TRUE)
  ec_counts <- ec_counts[o]
  ec_multimapping <- ec_multimapping[o]

  barplot(1 + ec_counts, col = ifelse(ec_counts > threshold * nrow(X), "blue", "red"), border = ifelse(ec_multimapping, par("bg"), "black"), lwd = 2, log = "y", las = 3)
  abline(h = 1 + threshold * nrow(X), lty = 2)

  return(names(ec_counts[ec_counts > threshold * nrow(X) & !ec_multimapping]))
}
ecidx <- plot_gene_filt("PTPRC")


## =============================================================================
## 2. Fitting multivariate logistic regression
## =============================================================================

## -----------------------------------------------------------------------------
## Load cell types
labels <- readLines("./cell.labels")

## Subset data matrices to match each cell type
X_naive <- X[labels == "Naive", ]
X_mem <- X[labels == "Mem", ]
X_cyto <- X[labels == "Cyto", ]

X_tx_naive <- X_tx[labels == "Naive", cd45_txs]
X_tx_mem <- X_tx[labels == "Mem", cd45_txs]
X_tx_cyto <- X_tx[labels == "Cyto", cd45_txs]


## -----------------------------------------------------------------------------
## Define wrappers to perform mLR and MAST.
## Adapted from code block 37 of:
## https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb

## For a set of TCCs (ecidx) and two data matrices (X1 and X2),
## perform multivariate logistic regression of X1[ecidx,] vs X2[ecidx,]
logr_ecidx <- function(ecidx, X1, X2) {
  library(lmtest)
  logr_labels <- rep(1:0, times = c(nrow(X1), nrow(X2))) ## Define dependent variable for the log.reg. models
  data <- as.data.frame(as.matrix(Matrix(rbind(X1[, ecidx], X2[, ecidx]), sparse = FALSE)))
  colnames(data) <- paste0("X", colnames(data)) ## Avoids issues in R in formulas with numeric column names
  data <- cbind(y = logr_labels, data)
  fmla <- paste0("y~", paste(colnames(data)[-1], collapse = "+"))
  glm <- glm(fmla, data = data, family = binomial()) ## Log.reg. model
  p <- lrtest(glm, glm(logr_labels ~ 1, family = binomial()))[2, 5] ## Test log.reg. vs null model
  p ## Return p-value
}

## Same as above, but using MAST instead of multivariate LR
mast_ecidx <- function(ecidx, X1, X2) {
  library(MAST)
  d <- t(as.matrix(rbind(X1[, ecidx], X2[, ecidx]))) ## Prepare data matrix for MAST
  rownames(d) <- paste0("X", rownames(d)) ## Makre sure feature names are not numeric, which would cause problems in R formulas

  gene_name <- unlist(EC_to_gene_names[ecidx]) ## Specify gene names
  if (length(gene_name) == 0) {
    gene_name <- rep("CD45", length(ecidx)) ## With default to CD45
  }

  ## Transform data matrix to single cell assay
  sca <- FromMatrix(
    log2(1 + d),
    cData = data.frame(row.names = colnames(d), y = rep(1:0, times = c(nrow(X1), nrow(X2))), cnecon = c(Matrix:::rowSums(X1 == 0), Matrix:::rowSums(X2 == 0))), ## cnecon: equivalent of the Cellular Detection Rate (Finak et al, 2015) at the TCC level
    fData = data.frame(row.names = paste0("X", ecidx), name = gene_name)
  )
  fmla <- "~y+scale(cnecon)" ## Model is TCC ~ cell type + CDR

  LRT <- "y"
  zlmCond <- zlm(as.formula(fmla), sca) ## Fit MAST
  summaryCond <- summary(zlmCond, doLRT = LRT) ## Get p-values for cell-type per TCC
  summaryDt <- summaryCond$datatable

  ## Reshape results: keep Hurdle p-values
  fcHurdle <- merge(
    summaryDt[contrast == LRT & component == "H", .(primerid, `Pr(>Chisq)`)], ## hurdle P values
    summaryDt[contrast == LRT & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
    by = "primerid" ## logFC coefficients
  )

  ## Correct
  fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, method = "fdr")
  ## Append TCC to genes mapping
  fcHurdle <- left_join(fcHurdle, as_tibble(rowData(sca)), by = c("primerid" = "primerid"))

  tx <- setNames(fcHurdle$`Pr(>Chisq)`, fcHurdle$primerid) ## TCC level p-values

  ## Split p-values (split = genes)
  fcHurdle <- split(fcHurdle, fcHurdle$name)
  ## Perform Sidak aggregation
  fcHurdle2 <- lapply(fcHurdle, function(slice) {
    n <- nrow(slice)
    x <- slice[, "Pr(>Chisq)"]
    x[x == 0] <- .Machine$double.xmin
    res <- 1 - (1 - min(x))^n ## First-order approximation for very low values
    if (res == 0) {
      res <- n * min(x)
    }
    res <- data.frame(name = slice[1, "name"], sidak = res, n = n, stringsAsFactors = FALSE)
  })
  fcHurdle2 <- do.call(rbind, fcHurdle2[order(names(fcHurdle2))])

  return(list(gn = setNames(fcHurdle2$sidak, fcHurdle2$primerid), tx = tx, fcHurdle = fcHurdle))
}


## -----------------------------------------------------------------------------
## This code block performs 200 comparisons for each pair of celltype.
## For each comparison, a random subsample of size 2,000 (out of 3,000)
## is drawn per cell type and compaired to a similar subsample for the
## other cell type, then we return the p-value for the four methods shown
## in that figure (univariate LR with Bonferroni correction, LR at the gene level,
## multivariate LR, and MAST-Sidak).
## We also investigate transcript-level analysis using mLR and MAST-Sidak
## (but the results contradict biology for both methods...).
mats <- list(X_mem = X_mem[, ecidx, drop = FALSE],
             X_naive = X_naive[, ecidx, drop = FALSE],
             X_cyto = X_cyto[, ecidx, drop = FALSE])

mats_tx <- list(X_tx_mem = X_tx_mem[, cd45_txs, drop = FALSE],
                X_tx_naive = X_tx_naive[, cd45_txs, drop = FALSE],
                X_tx_cyto = X_tx_cyto[, cd45_txs, drop = FALSE])

file <- "./rdata/Fig2.RData"
if (!file.exists(file)) {
  set.seed(123)
  ps <- lapply(
    combn(1:3, 2, simplify = FALSE),
    function(indices) {
      res <- sapply(
        1:200,
        function(x) {
          print(x)
          n <- 2000
          s1 <- sample(1:3000, n)
          s2 <- sample(1:3000, n)
          y <- rep(1:0, each = n)

          ## Log. reg. on TCC
          logreg_tcc <- logr_ecidx(ecidx, mats[[indices[1]]][s1, ecidx], mats[[indices[2]]][s2, ecidx])

          ## MAST on TCC
          mast_tcc <- mast_ecidx(ecidx, mats[[indices[1]]][s1, ecidx], mats[[indices[2]]][s2, ecidx])$gn 

          ## TCC to genes aggregation
          tcc <- as.matrix(t(rbind(mats[[indices[1]]][s1, ecidx], mats[[indices[2]]][s2, ecidx])))
          rownames(tcc) <- paste0("X", rownames(tcc))
          gn <- colSums(tcc)

          ## Log. reg. on tcc "Independent tests", Gene counts, transcripts (tx), and mast on transcripts
          logreg_indep <- min(p.adjust(
              do_logreg(tcc, y, tx2gene = rownames(tcc), covariates = NULL),
              n = length(ecidx), method = "bonferroni")) 
          logreg_gene <- lrtest(glm(y ~ gn, family = binomial()),
                                glm(y ~ 1, family = binomial()))[2, 5]
          logreg_tx <- logr_ecidx(cd45_txs,
                                  mats_tx[[indices[1]]][s1, cd45_txs],
                                  mats_tx[[indices[2]]][s2, cd45_txs])
          mast_tx <- mast_ecidx(cd45_txs,
                                mats_tx[[indices[1]]][s1, cd45_txs],
                                mats_tx[[indices[2]]][s2, cd45_txs])$gn

          return(c(logreg_tcc, mast_tcc, logreg_indep, logreg_gene, logreg_tx, mast_tx))
        }
      )
      rownames(res) <- c("logreg_tcc", "mast_tcc", "logreg_indep", "logreg_gene", "logreg_tx", "mast_tx")
      res
    }
  )
  names(ps) <- sapply(
    combn(1:3, 2, simplify = FALSE),
    function(x) {
      paste0(names(mats[x]), collapse = " vs ")
    }
  )
  save(ps, file = file)
} else {
  load(file)
}



## =============================================================================
## 3. Plots
## =============================================================================

## -----------------------------------------------------------------------------
## Includes all methods mentioned prior (figure not included in paper)
## Plot as histogram with overlaid density and vertical line at the median
png("Fig 2_EB.png", width = 1000, height = 3000, res = 300)
par(mfrow = c(3, 1), mar = c(3, 3, 3, 1), cex.axis = 1.2, cex.lab = 1.2)
lapply(
  names(ps),
  function(i) {
    elt <- ps[[i]]
    elt <- elt[, apply(elt, 2, min) > 0]

    fig2a <- as.list(as.data.frame(t(elt)))
    fig2a <- lapply(fig2a, function(x) {
      x <- -log10(x)
      x
    })
    xlim <- range(sapply(fig2a, range))
    densities <- lapply(fig2a, density, from = 0, to = ceiling(xlim[2]), n = 1024)
    hists <- sapply(
      names(fig2a),
      function(x) {
        hist(fig2a[[x]], plot = FALSE)
      },
      simplify = FALSE
    )
    ylim <- range(unlist(sapply(densities, function(x) {
      x$y
    })), unlist(sapply(hists, function(x) {
      x$density
    })))
    ylim[2] <- ceiling(ylim[2] * 10 + 1) / 10

    library(gplots)
    cols <- sapply(c("logreg_gene" = "red",
                     "logreg_indep" = "blue",
                     "logreg_tcc" = "orange",
                     "mast_tcc" = "purple",
                     "logreg_tx" = "green",
                     "mast_tx" = "darkgreen"),
                   col2hex)
    tp <- function(hexcol, tp = "44") {
      setNames(paste0(substr(hexcol, 1, 7), tp), names(hexcol))
    }
    cols <- setNames(tp(cols, "99"), names(cols))

    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    box(bty = "l")
    sapply(
      names(fig2a),
      function(x) {
        plot(hists[[x]], border = cols[x], col = tp(cols[x]), add = TRUE, freq = FALSE)
      }
    )
    sapply(
      names(fig2a),
      function(x) {
        abline(v = median(fig2a[[x]]), lty = 2, lwd = 2, col = cols[x])
      }
    )
    sapply(
      names(fig2a),
      function(x) {
        lines(densities[[x]], col = tp(cols[x], "FF"), lwd = 2)
      }
    )
    axis(side = 1, at = seq(xlim[1], xlim[2], by = 1))
    title(xlab = "-log10(P)", main = tools::toTitleCase(gsub("X_", "", i)))
    axis(side = 2)

    legend(x = "topright", lty = 1, legend = names(cols), col = cols, bty = "n", cex = 1.2)
  }
)
dev.off()


## -----------------------------------------------------------------------------
## Simplified figure included in paper using select methods 
library(RSvgDevice)
## Select methods
ps <- lapply(ps, function(x) {
  x[c("logreg_tcc", "mast_tcc", "logreg_indep", "logreg_gene"), ]
})
## Define bin boundaries for histograms
breaks <- lapply(ps, range)
breaks <- lapply(breaks, function(x) {
  seq(-log10(x[1]), -log10(x[2]), length.out = 20)
})
svg("Fig CD45.svg", width = 3, height = 9)
par(mfrow = c(3, 1), mar = c(3, 3, 3, 1), cex.axis = 1.2, cex.lab = 1.2)
lapply(
  names(ps),
  function(i) {
    elt <- ps[[i]]

    ## -log10 transform
    fig2a <- as.list(as.data.frame(t(elt)))
    fig2a <- lapply(fig2a, function(x) {
      x <- -log10(x)
      x
    })

    ## Compute densities
    xlim <- range(sapply(fig2a, range))
    densities <- lapply(fig2a, density, from = 0, to = ceiling(xlim[2]), n = 1024)

    ## Compute histograms
    hists <- sapply(
      names(fig2a),
      function(x) {
        hist(fig2a[[x]], plot = FALSE, breaks = breaks[[i]])
      },
      simplify = FALSE
    )

    ## Find y-axis ranges
    ylim <- range(unlist(sapply(densities, function(x) {
      x$y
    })), unlist(sapply(hists, function(x) {
      x$density
    })))
    ylim[2] <- ceiling(ylim[2] * 10 + 1) / 10

    ## Choose colors and transparency
    library(gplots)
    cols <- sapply(c("logreg_gene" = "red", "logreg_indep" = "blue", "logreg_tcc" = "orange", "mast_tcc" = "purple"), col2hex)
    tp <- function(hexcol, tp = "44") {
      setNames(paste0(substr(hexcol, 1, 7), tp), names(hexcol))
    }
    cols <- setNames(tp(cols, "CC"), names(cols))

    ## Plot
    plot.new()
    plot.window(xlim = xlim, ylim = ylim)
    box(bty = "l")

    ## Histograms
    sapply(
      names(fig2a),
      function(x) {
        plot(hists[[x]], border = cols[x], col = tp(cols[x]), add = TRUE, freq = FALSE, lwd = 2)
      }
    )

    ## Median dashed line
    sapply(
      names(fig2a),
      function(x) {
        abline(v = median(fig2a[[x]]), lty = 2, lwd = 2, col = cols[x])
      }
    )

    ## Densities
    sapply(
      names(fig2a),
      function(x) {
        lines(densities[[x]], col = tp(cols[x], "FF"), lwd = 2)
      }
    )

    ## Add annotations
    axis(side = 1, at = seq(xlim[1], xlim[2], by = 1))
    title(xlab = "-log10(P)", main = tools::toTitleCase(gsub("X_", "", i)))
    axis(side = 2)
    legend(x = "topright", lty = 1, legend = names(cols), col = cols, bty = "n", cex = 1.2)
  }
)
dev.off()
