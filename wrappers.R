library(lmtest)
library(parallel)
library(MAST)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)

#################### This script contains wrappers used throughout other scripts, mostly to run
#################### various Gene Differential Expression methods

## Wrapper for logistic regression
do_logreg <- function(data, response, tx2gene, covariates = NULL) {
    stopifnot(nrow(data) == length(tx2gene))
    
    ## Split into a list of matrix. Each element = matrix of (transcripts mapping to a
    ## given gene x cells)
    y <- response
    tx_split <- split_matrix(as.matrix(data), tx2gene)
    tx_split <- tx_split[sort(unique(tx2gene))]
    tx_split <- lapply(tx_split, t)
    tx_split <- lapply(tx_split, cbind, y)
    if (!is.null(colnames(covariates))) {
        tx_split <- lapply(tx_split, cbind, covariates)  ## Adding potential covariates (e.g. batch effect, CDR)
    }
    tx_split <- lapply(tx_split, as.data.frame)
    
    ## On each gene-wise matrix, run multivariate LR and report likelihood-ratio test
    ## p-value
    regs <- mclapply(tx_split, function(x, y) {
        fmla <- paste("y~", paste(setdiff(colnames(x), "y"), collapse = "+"), sep = "")
        glm <- glm(fmla, data = x, family = binomial())
        if (is.null(colnames(covariates))) {
            nullmodel <- glm(x[, "y"] ~ 1, data = x, family = binomial())
        } else {
            fmla <- paste0("y~", paste0(colnames(covariates), collapse = "+"))
            nullmodel <- glm(fmla, data = x, family = binomial())
        }
        p <- lrtest(glm, nullmodel)[2, 5]
    }, y = y, mc.cores = getOption("mc.cores", 1L))
    
    return(setNames(unlist(regs), names(tx_split)))
}

## Wrapper for MAST - Sidak data: sca=SingleCellAssay with all the annotations
do_mast_tx <- function(fmla, sca, LRT) {
    
    ## Make sure we are splitting the matrix using a character vector. Factors would
    ## lead in extremely slow behaviour.
    rowData(sca)[, "gene_ensembl"] <- as.character(rowData(sca)[, "gene_ensembl"])
    
    ## Run MAST
    zlmCond <- zlm(as.formula(fmla), sca)
    summaryCond <- summary(zlmCond, doLRT = LRT)  ## Get p-values for cell types at the transcript level
    summaryDt <- summaryCond$datatable
    
    ## Extract results and reshape
    fcHurdle <- merge(summaryDt[contrast == LRT & component == "H", .(primerid, `Pr(>Chisq)`)], 
        summaryDt[contrast == LRT & component == "logFC", .(primerid, coef, ci.hi, 
            ci.lo)], by = "primerid"  # logFC coefficients
)
    summaryDt_LRT <- subset(summaryDt, contrast == LRT & component %in% c("H", "D", 
        "C"))
    summaryDt_LRT <- left_join(summaryDt_LRT, as_tibble(rowData(sca)), by = c(primerid = "transcript"))
    summaryDt_LRT <- split(summaryDt_LRT, summaryDt_LRT$component)
    
    ## Sidak aggregation of p-values from transcripts to genes
    summaryDt_LRT <- lapply(summaryDt_LRT, function(x) {
        do.call(rbind, mclapply(split(x, as.character(x$gene_ensembl)), FUN = function(slice) {
            n <- nrow(slice)
            x <- slice[, "Pr(>Chisq)"]
            x[x == 0] <- .Machine$double.xmin
            res <- 1 - (1 - min(x))^n
            if (res == 0) {
                res <- n * min(x)
            }
            res <- data.frame(gene_ensembl = slice[1, "gene_ensembl"], sidak = res, 
                n = n, stringsAsFactors = FALSE)
            res
        }, mc.cores = getOption("mc.cores", 1L)))
    })
    
    fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, method = "fdr")
    fcHurdle <- left_join(fcHurdle, as_tibble(rowData(sca)), by = c(primerid = "transcript"))
    
    res <- list()
    res$tx <- setNames(fcHurdle$`Pr(>Chisq)`, fcHurdle$primerid)  ## Save p-values MAST univariate on transcripts
    res$fcHurdle <- fcHurdle
    res$summaryDt <- summaryDt
    res$summaryDt_LRT <- summaryDt_LRT
    res$gn_list <- lapply(summaryDt_LRT, function(x) {
        setNames(x[, "sidak"], x[, "gene_ensembl"])
    })
    
    ## This should give identical result as the code that is directly below yet it
    ## outputs many NaN (Inf) when there is no apparent reason to. So we rewrote this
    ## in base R instead of tidyverse
    ## fcHurdle%>%as_tibble()%>%group_by(gene_ensembl)%>%summarise( sidak={
    ## x=get('Pr(>Chisq)') x[x==0]=.Machine$double.xmin res=1-(1-min(x))^n()
    ## if(res==0){ res=n*min(x) } res }, n=n() )->fcHurdle2
    
    ## Perform MAST univariate on genes
    fcHurdle <- split(fcHurdle, as.character(fcHurdle$gene_ensembl))
    fcHurdle2 <- lapply(fcHurdle, function(slice) {
        n <- nrow(slice)
        x <- slice[, "Pr(>Chisq)"]
        x[x == 0] <- .Machine$double.xmin
        res <- 1 - (1 - min(x))^n
        if (res == 0) {
            res <- n * min(x)
        }
        res <- data.frame(gene_ensembl = slice[1, "gene_ensembl"], sidak = res, n = n, 
            stringsAsFactors = FALSE)
    })
    fcHurdle2 <- do.call(rbind, fcHurdle2[order(names(fcHurdle2))])
    
    res$gn <- setNames(fcHurdle2$sidak, fcHurdle2$gene_ensembl)
    return(res)
}

## Not used. Would do MAST univariate on genes. We use do_mast_tx() instead and
## use the $tx element of the output to obtain the same result (although it is
## slower) do_mast_gn=function(fmla,sca,LRT){ ##data: SingleCellAssay with all the
## annotations zlmCond=zlm(as.formula(fmla),sca) summaryCond=summary(zlmCond,
## doLRT=LRT) summaryDt=summaryCond$datatable

## fcHurdle = merge( summaryDt[contrast==LRT & component=='H',.(primerid,
## `Pr(>Chisq)`)], #hurdle P values summaryDt[contrast==LRT & component=='logFC',
## .(primerid, coef, ci.hi, ci.lo)], by='primerid' #logFC coefficients )

## fcHurdle$fdr = p.adjust(fcHurdle$`Pr(>Chisq)`,method='fdr')

## return(list(gn=setNames(fcHurdle$`Pr(>Chisq)`,fcHurdle$primerid),fcHurdle=fcHurdle))
## }

## Sidak aggregation of p-values
sidak <- function(x) {
    n <- length(x)
    x[x == 0] <- .Machine$double.xmin
    res <- 1 - (1 - min(x))^n
    if (res == 0) {
        res <- n * min(x)
    }
    res
}

do_t_tests <- function(response, data, t2g, covariates = NULL) {
    if (missing(covariates)) {
        covariates <- rep(1, length(response))
    }
    response <- as.factor(response)
    aovs <- apply(data, 1, function(x, response, covariates) {
        summary(aov(x ~ response + covariates))[[1]][[5]][1]
    }, response = response, covariates = covariates)
    res <- list()
    res$tx <- aovs
    res$tx[is.nan(res$tx)] <- 1
    res$gn <- sapply(split(res$tx, t2g), sidak)
    res
}

plot_gene <- function(sca) {
    gene_data <- as.data.table(rowData(sca))
    cell_data <- as.data.table(colData(sca))
    gene_data[, `:=`(transcript, paste0(transcript, ifelse(de_transcript, " (DE)", 
        "")))]
    rownames(sca) <- unlist(gene_data[, "transcript"])
    
    dt <- data.table(assays(sca)$tpm, keep.rownames = "transcript")
    dt_long <- suppressWarnings(melt(dt, variable.name = "Cell", value.name = "tpm"))
    dt_long <- merge(gene_data, dt_long, by = "transcript")
    dt_long <- merge(cell_data, dt_long, by = "Cell")
    
    ggplot(dt_long, aes(x = paste0(sub("Batch", "B", Batch), "_", sub("Group", "G", 
        Group)), y = log10(1 + tpm))) + geom_violin() + geom_jitter() + facet_wrap(~transcript) + 
        theme_bw() + xlab(paste0(unique(gene_data$gene_ensembl), collapse = "/"))
}

## Adapted from
## https://github.com/pachterlab/NYMP_2018/blob/master/simulations/RSEM/R/roc_helpers.R
## should be slightly faster and also returns false positive rate
library(dplyr)
calculate_fdr <- function(true_DE, pvalues, title) {
    df <- data.frame(true_DE = true_DE, pvalues = pvalues)
    df <- df[order(df$pvalues), ]
    total_positive <- sum(true_DE)
    total_negative <- sum(!true_DE)
    df <- df %>% group_by(pvalues) %>% summarise(p = sum(true_DE), n = sum(!true_DE), 
        l = n())
    ## print(head(df)) fdr <- sapply(seq(df$pvalues), function(i) sum(df$n[1:i]) /
    ## sum(df$l[1:i]))
    fdr <- cumsum(df$n)/cumsum(df$l)
    ## sensitivity <- sapply(seq(df$pvalues), function(i)
    ## sum(df$p[1:i])/total_positive)
    sensitivity <- cumsum(df$p)/total_positive
    ## falsepositiverate <- sapply(seq(df$pvalues), function(i)
    ## sum(df$n[1:i])/total_negative)
    falsepositiverate <- cumsum(df$n)/total_negative
    
    n_found <- sapply(seq(df$pvalues), function(i) sum(df$p[1:i]))
    n_de <- sapply(seq(df$pvalues), function(i) sum(df$l[1:i]))
    five <- min(which(df$pvalues > 0.05)) - 1
    ten <- min(which(df$pvalues > 0.1)) - 1
    twenty <- min(which(df$pvalues > 0.2)) - 1
    
    list(fdr = fdr, sensitivity = sensitivity, n_found = n_found, n_de = n_de, pvalues = df$pvalues, 
        five = five, ten = ten, twenty = twenty, title = title, FalsePositiveRate = falsepositiverate)
}


#### Below are misc. functions used throughout the other scripts:

#' Aggregates matrix or data frame
#'
#' @param df A matrix or data frame
#' @param groups A vector of groups (discrete values)
#' @param fun A function that aggregates a vector into objects of length 1
#' @param margin If 1, aggregates rows, if 2 aggregates columns. Defaults to 1
#' @param ... passed to fun
#' @return A data.frame with aggregated rows or columns
#' @export

aggregate_df <- function(df, groups, fun = mean, margin = 1, ...) {
    if (length(groups) != dim(df)[margin]) {
        stop("Size of 'groups' vector is different that the one of the specified data margin")
    }
    if (is.data.frame(df)) {
        if (margin == 2) {
            df <- as.data.frame(t(df))
        } else {
            df <- as.data.frame(df)
        }
        df <- split(df, groups)
    } else if (is.matrix(df)) {
        df <- split_matrix(df, groups, byrow = margin == 1)
    }
    
    res <- do.call(rbind, lapply(df, function(x) {
        apply(x, 2, fun, ...)
    }))
    
    if (margin == 2) {
        return(t(res))
    } else {
        return(res)
    }
}

#' Makes a color transparent
#' @export
#' @param colors A vector of colors as in `?col2rgb`
#' @param alpha transparency value (0=fully transparent, 255=fully opaque)
#'
## credit
## http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color

makeTransparent <- function(colors, alpha = 255) {
    sapply(colors, function(col) {
        col <- col2rgb(col)
        rgb(red = col[1, ], green = col[2, ], blue = col[3, ], alpha = alpha, maxColorValue = 255)
    })
}

#' Fast splitting of matrix to list (avoids conversion to data.frame)
#' @export
split_matrix <- function(mat, vector, byrow = TRUE) {
    if (byrow & nrow(mat) != length(vector)) {
        stop("if byrow=TRUE, vector's length should have length nrow(mat)")
    } else if (!byrow & ncol(mat) != length(vector)) {
        !byrow & ncol(mat) != length(vector)
        stop("if byrow=FALSE, vector's length should have length ncol(mat)")
    }
    
    if (byrow) {
        levels <- split(1:nrow(mat), vector)
        res <- lapply(levels, function(x) {
            mat[x, , drop = FALSE]
        })
    } else {
        levels <- split(1:ncol(mat), vector)
        res <- lapply(levels, function(x) {
            mat[, x, drop = FALSE]
        })
    }
    res
}
