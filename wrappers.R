library(lmtest)
library(parallel)
library(ebmisc)
library(MAST)
library(MultiAssayExperiment)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)

####################
## This script contains wrappers used throughout other scripts, mostly to run various Gene Differential Expression methods
####################

## Wrapper for logistic regression
do_logreg=function(data,response,tx2gene,covariates=NULL){
    stopifnot(nrow(data)==length(tx2gene))

    ## Split into a list of matrix. Each element = matrix of (transcripts mapping to a given gene x cells)
    y = response
    if(require("ebmisc")){
        tx_split=split_matrix(as.matrix(data),tx2gene)
    } else {
        tx_split=split(as.data.frame(as.matrix(data)),tx2gene)
    }
    tx_split=tx_split[sort(unique(tx2gene))]
    tx_split=lapply(tx_split,t)
    tx_split=lapply(tx_split,cbind,y)
    if(!is.null(colnames(covariates))){
        tx_split=lapply(tx_split,cbind,covariates) ## Adding potential covariates (e.g. batch effect, CDR)
    }
    tx_split=lapply(tx_split,as.data.frame)

    ## On each gene-wise matrix, run multivariate LR and report likelihood-ratio test p-value
    regs=mclapply(
        tx_split,
        function(x,y){
            fmla=paste("y~",paste(setdiff(colnames(x),"y"),collapse="+"),sep="")
            glm=glm(fmla,data=x,family=binomial())
            if(is.null(colnames(covariates))){
                nullmodel=glm(x[,"y"]~1,data=x,family=binomial())
            } else {
                fmla=paste0("y~",paste0(colnames(covariates),collapse="+"))
                nullmodel=glm(fmla,data=x,family=binomial())
            }
            p=lrtest(glm,nullmodel)[2,5]
        },
        y=y,
        mc.cores=getOption("mc.cores",1L)
    )
    
    return(setNames(unlist(regs),names(tx_split)))
}

## Wrapper for MAST - Sidak
do_mast_tx=function(fmla,sca,LRT){ ##data: SingleCellAssay with all the annotations

    ## Make sure we are splitting the matrix using a character vector. Factors would lead in extremely slow behaviour.
    rowData(sca)[,"gene_ensembl"]=as.character(rowData(sca)[,"gene_ensembl"])

    ## Run MAST
    zlmCond=zlm(as.formula(fmla),sca)
    summaryCond=summary(zlmCond, doLRT=LRT) ## Get p-values for cell types at the transcript level
    summaryDt=summaryCond$datatable

    ## Extract results and reshape
    fcHurdle = merge(
        summaryDt[contrast==LRT & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
        summaryDt[contrast==LRT & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid' #logFC coefficients
    )
    summaryDt_LRT=subset(summaryDt,contrast==LRT&component%in%c("H","D","C"))
    summaryDt_LRT=left_join(summaryDt_LRT,as_tibble(rowData(sca)),by=c("primerid"="transcript"))
    summaryDt_LRT=split(summaryDt_LRT,summaryDt_LRT$component)

    ## Sidak aggregation of p-values from transcripts to genes
    summaryDt_LRT=lapply(
        summaryDt_LRT,
        function(x){
            do.call(
                rbind,
                mclapply(
                    split(x,as.character(x$gene_ensembl)),
                    FUN=function(slice){
                        n=nrow(slice)
                        x=slice[,"Pr(>Chisq)"]
                        x[x==0]=.Machine$double.xmin
                        res=1-(1-min(x))^n
                        if(res==0){
                            res=n*min(x)
                        }
                        res=data.frame(gene_ensembl=slice[1,"gene_ensembl"],sidak=res,n=n,stringsAsFactors=FALSE)
                        res
                    },
                    mc.cores=getOption("mc.cores",1L)
                )
            )
        }
    )
    
    fcHurdle$fdr = p.adjust(fcHurdle$`Pr(>Chisq)`,method="fdr")
    fcHurdle=left_join(fcHurdle,as_tibble(rowData(sca)),by=c("primerid"="transcript"))
    
    res=list()
    res$tx=setNames(fcHurdle$`Pr(>Chisq)`,fcHurdle$primerid) ## Save p-values MAST univariate on transcripts
    res$fcHurdle=fcHurdle
    res$summaryDt=summaryDt
    res$summaryDt_LRT=summaryDt_LRT
    res$gn_list=lapply(summaryDt_LRT,function(x){
        setNames(x[,"sidak"],x[,"gene_ensembl"])
    })
    
    ## This should give identical result as the code that is directly below yet it outputs many NaN (Inf) when there is no apparent reason to. So we rewrote this in base R instead of tidyverse
    ## fcHurdle%>%as_tibble()%>%group_by(gene_ensembl)%>%summarise(
    ##                                                       sidak={
    ##                                                           x=get("Pr(>Chisq)")
    ##                                                           x[x==0]=.Machine$double.xmin
    ##                                                           res=1-(1-min(x))^n()
    ##                                                           if(res==0){
    ##                                                               res=n*min(x)
    ##                                                           }
    ##                                                           res
    ##                                                       },
    ##                                                       n=n()
    ##                                                   )->fcHurdle2

    ## Perform MAST univariate on genes
    fcHurdle=split(fcHurdle,as.character(fcHurdle$gene_ensembl))
    fcHurdle2=lapply(fcHurdle,function(slice){
        n=nrow(slice)
        x=slice[,"Pr(>Chisq)"]
        x[x==0]=.Machine$double.xmin
        res=1-(1-min(x))^n
        if(res==0){
            res=n*min(x)
        }
        res=data.frame(gene_ensembl=slice[1,"gene_ensembl"],sidak=res,n=n,stringsAsFactors=FALSE)
    })
    fcHurdle2=do.call(rbind,fcHurdle2[order(names(fcHurdle2))])

    res$gn=setNames(fcHurdle2$sidak,fcHurdle2$gene_ensembl)
    return(res)
}

## Not used. Would do MAST univariate on genes. We use do_mast_tx() instead and use the $tx element of the output to obtain the same result (although it is slower)
## do_mast_gn=function(fmla,sca,LRT){ ##data: SingleCellAssay with all the annotations
##     zlmCond=zlm(as.formula(fmla),sca)
##     summaryCond=summary(zlmCond, doLRT=LRT)
##     summaryDt=summaryCond$datatable

##     fcHurdle = merge(
##         summaryDt[contrast==LRT & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
##         summaryDt[contrast==LRT & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid' #logFC coefficients
##     )

##     fcHurdle$fdr = p.adjust(fcHurdle$`Pr(>Chisq)`,method="fdr")
    
##     return(list(gn=setNames(fcHurdle$`Pr(>Chisq)`,fcHurdle$primerid),fcHurdle=fcHurdle))
## }

## Sidak aggregation of p-values
sidak=function(x){
    n=length(x)
    x[x==0]=.Machine$double.xmin
    res=1-(1-min(x))^n
    if(res==0){
        res=n*min(x)
    }
    res
}

do_t_tests=function(response,data,t2g,covariates=NULL){
    if(missing(covariates)){
        covariates=rep(1,length(response))
    }
    response=as.factor(response)
    aovs=apply(
        data,
        1,
        function(x,response,covariates){
            summary(aov(x~response+covariates))[[1]][[5]][1]
        },
        response=response,covariates=covariates
    )
    res=list()
    res$tx=aovs
    res$tx[is.nan(res$tx)]=1
    res$gn=sapply(split(res$tx,t2g),sidak)
    res
}

plot_gene=function(sca){
    gene_data=as.data.table(rowData(sca))
    cell_data=as.data.table(colData(sca))
    gene_data[,transcript:=paste0(transcript,ifelse(de_transcript," (DE)",""))]
    rownames(sca)=unlist(gene_data[,"transcript"])

    dt=data.table(assays(sca)$tpm,keep.rownames="transcript")
    dt_long=suppressWarnings(melt(dt,variable.name="Cell",value.name="tpm"))
    dt_long=merge(gene_data,dt_long,by="transcript")
    dt_long=merge(cell_data,dt_long,by="Cell")

    ggplot(dt_long, aes(x = paste0(sub("Batch","B",Batch),"_",sub("Group","G",Group)), y = log10(1+tpm))) + 
        geom_violin() + 
        geom_jitter() +
        facet_wrap(~transcript) +
        theme_bw() +
        xlab(paste0(unique(gene_data$gene_ensembl),collapse="/"))
}

## Adapted from https://github.com/pachterlab/NYMP_2018/blob/master/simulations/RSEM/R/roc_helpers.R should be slightly faster and also returns false positive rate
library(dplyr)
calculate_fdr <- function(true_DE, pvalues, title)
{
    df <- data.frame(true_DE = true_DE, pvalues=pvalues)
    df <- df[order(df$pvalues), ]
    total_positive <- sum(true_DE)
    total_negative <- sum(!true_DE)
    df <- df %>% group_by(pvalues) %>% summarise(p = sum(true_DE), n = sum(!true_DE), l=n())
    ## print(head(df))
    ## fdr <- sapply(seq(df$pvalues), function(i) sum(df$n[1:i]) / sum(df$l[1:i]))
    fdr = cumsum(df$n)/cumsum(df$l)
    ## sensitivity <- sapply(seq(df$pvalues), function(i) sum(df$p[1:i])/total_positive)
    sensitivity = cumsum(df$p)/total_positive
    ## falsepositiverate <- sapply(seq(df$pvalues), function(i) sum(df$n[1:i])/total_negative)
    falsepositiverate = cumsum(df$n)/total_negative
    
    n_found <- sapply(seq(df$pvalues), function(i) sum(df$p[1:i]))
    n_de <- sapply(seq(df$pvalues), function(i) sum(df$l[1:i]))
    five <- min(which(df$pvalues > .05)) -1 
    ten <- min(which(df$pvalues > .1)) - 1
    twenty <- min(which(df$pvalues > .2)) -1
    
    list(
        fdr = fdr,
        sensitivity = sensitivity,
        n_found = n_found,
        n_de = n_de,
        pvalues = df$pvalues,
        five=five, ten=ten, twenty=twenty, title=title,
        FalsePositiveRate=falsepositiverate
    )
}

sidak = function(slice){
    n=length(slice)
    x=slice
    x[is.na(x)]=1
    x[x==0]=.Machine$double.xmin
    res=1-(1-min(x))^n
    if(res==0){
        res=n*min(x)
    }
    res
}

##Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
do_edgeRCDR = function(counts, groups){
    require(edgeR)
    cdr <- scale(colMeans(counts > 0))
    design <- model.matrix(~ cdr + groups)
    dge <- DGEList(counts, group = groups)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)
    tt = tt@.Data[[1]]
    edgeRQLF_CDR = setNames(tt[, "PValue"], rownames(tt))
    edgeRQLF_CDR
}

##Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
do_edgeR = function(counts, groups){
    require(edgeR)
    design <- model.matrix(~ groups)
    dge <- DGEList(counts, group = groups)
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design = design)
    fit <- glmQLFit(dge, design = design)
    qlf <- glmQLFTest(fit)
    tt <- topTags(qlf, n = Inf)
    tt = tt@.Data[[1]]
    edgeRQLF = setNames(tt[, "PValue"], rownames(tt))
    edgeRQLF
}

##Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
do_limma = function(counts, groups){
    require(limma)
    require(edgeR)
    ## limmatrend
    design <- model.matrix(~ groups)
    dge <- DGEList(counts, group = groups)
    dge <- calcNormFactors(dge)
    y <- new("EList")
    y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
    fit <- lmFit(y, design = design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    tt <- topTable(fit, n = Inf, adjust.method = "BH")
    limmatrend = setNames(tt[, "P.Value"], rownames(tt))
    limmatrend
}

##Adapted from https://github.com/csoneson/conquer_comparison/tree/master/scripts
do_DESeq2 = function(counts, groups){
    require(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = round(counts + 1), 
                                  colData = data.frame(condition = groups), 
                                  design = ~condition)
    dds <- DESeq(dds, betaPrior = FALSE, fitType = "mean")
    ## res <- results(
    ##     dds,
    ##     contrast = c(
    ##         "condition",
    ##         levels(factor(groups))[1], 
    ##         levels(factor(groups))[2]),
    ##     alpha = 0.05
    ## )
    DESeq2 = setNames(
        rowData(dds)[,grep("WaldPvalue_condition", colnames(rowData(dds)), value = TRUE)],
        rownames(dds)
    )
    DESeq2
}

auc = function(TPR, FPR){
    stopifnot(length(TPR)==length(FPR))
    i = 2:length(TPR) - 1
    
    heights = (TPR[i] + TPR[i+1])/2
    widths = FPR[i+1]-FPR[i]
    sum(heights*widths)
}
