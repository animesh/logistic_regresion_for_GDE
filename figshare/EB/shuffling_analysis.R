########################################
## README:
##########
## This script
## 1. loads and harmize three real scRNAseq datasets for which multiple cell types are pre-defined
## 2. then perform GDE on each dataset (with two distinct expression frequency threshold for genes)
## 3. then shuffle cell type labels and perform GDE on shuffled labels (performed 5 times)
## 4. plot results
########################################

library(ggpubr)
library(splatter)
library(ebmisc)
library(MultiAssayExperiment)
library(Matrix)
library(hypergate)
library(rio)
library(data.table)
library(car)
library(RColorBrewer)
library(abind)

options("mc.cores"=1L,scipen=999) ## I have had issues running this script in parallel (memory issues, forking crashes) so I advise against mc.cores>1. scipen=999 is used to avoid exponential notation when converting numerics with as.character()

setwd("./figshare/EB/")
source("./wrappers.R")
splatter_dir="./figshare/splatter"

####################
## 1. Load and harmonize input.
####################

##
## 10X data (Zheng et al, Nat Com, 2017). Processed by Nnatros et al. Parsing adapted from https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
##

## Load data
X=readMM("./NYMP_2018/10x_example-logR/matrix.tcc.mtx")
colnames(X)=1:ncol(X)-1

## Mapping of TCCs to transcripts (EC_dict), TX numbers to ensembl gene transcript identifiers (TX_TO_ENST), to ensembl gene identifiers (TX_TO_ENSG)
TX_TO_ENST=read.csv("./NYMP_2018/10x_example-logR/TX_to_ENST.csv",sep="\t",stringsAsFactors=FALSE,header=FALSE)
TX_TO_ENST=setNames(TX_TO_ENST[,2],TX_TO_ENST[,1])
TX_TO_ENSG=read.csv("./NYMP_2018/10x_example-logR/TX_to_ENSG.csv",sep="\t",stringsAsFactors=FALSE,header=FALSE)
TX_TO_ENSG=setNames(TX_TO_ENSG[,2],TX_TO_ENSG[,1])
ENSG_to_name=read.csv("./NYMP_2018/10x_example-logR/ENSG_to_name.csv",sep="\t",stringsAsFactors=FALSE,header=FALSE)
ENSG_to_name=setNames(ENSG_to_name[,2],ENSG_to_name[,1])
EC_dict=read.csv("./NYMP_2018/10x_example-logR/matrix.ec",sep="\t",stringsAsFactors=FALSE,header=FALSE)
EC_dict=setNames(strsplit(EC_dict[,2],","),EC_dict[,1])
## Mapping TCCs to gene names and vice versa
EC_to_gene_names=stack(EC_dict,stringsAsFactors=FALSE)
colnames(EC_to_gene_names)=c("TX","TCC")
EC_to_gene_names=apply(EC_to_gene_names,2,as.character)
EC_to_gene_names=cbind(EC_to_gene_names,ENSG=TX_TO_ENSG[EC_to_gene_names[,"TX"]])
EC_to_gene_names=cbind(EC_to_gene_names,name=ENSG_to_name[EC_to_gene_names[,"ENSG"]])
gene_names_to_ECs=split(EC_to_gene_names[,"TCC"],EC_to_gene_names[,"name"])
gene_names_to_ECs=lapply(gene_names_to_ECs,unique)
gene_names_to_ECs=lapply(gene_names_to_ECs,sort)
EC_to_gene_names=split(EC_to_gene_names[,"name"],EC_to_gene_names[,"TCC"])
EC_to_gene_names=lapply(EC_to_gene_names,unique)
EC_to_gene_names=lapply(EC_to_gene_names,sort)

## Making sure mapping and matrix are in the same order. This is where scipen=999 is important.
tccs=intersect(names(EC_to_gene_names[sapply(EC_to_gene_names,length)==1]),colnames(X))
X=X[,tccs]
EC_to_gene_names=unlist(EC_to_gene_names[tccs])

##
## EMTAB3929 (Petropoulos et al)
##
EMTAB3929=readRDS("./data/EMTAB3929.rds") ## http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/EMTAB3929.rds
EMTAB3929=EMTAB3929[,grepl("E3",colnames(EMTAB3929)[["tx"]])|grepl("E4",colnames(EMTAB3929)$tx)] ## We compare E3 vs E4

##
## GSE64016
##
GSE64016=readRDS("./data/GSE64016.rds") ## http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/GSE64016.rds
GSE64016=GSE64016[,colData(GSE64016)$source_name_ch1%in%c("single H1-Fucci cell sorted from G1 phase of the cell cycle only","single H1-Fucci cell sorted from S phase of the cell cycle only")] ## We compare cells at the G1 vs S phases of the cell cycle

## List datasets and harmonize. List of list. Sublists elements are 1) xp : expression matrix, 2) grps : vector of cell types, 3) name : dataset's name, 4) t2g: Transcript (or TCCs) to gene mapping
data=list(
    list(
        xp=cbind(t(X[1:3000,]),t(X[3001:6000,])), ## Select CD4 Naive vs CD4 Memory only
        grps=rep(c(TRUE,FALSE),each=3000),
        name="Authors' 10X data",
        t2g=EC_to_gene_names
    ),
    list(
        xp=assays(EMTAB3929)$tx,
        grps=grepl("E3",colnames(EMTAB3929)[["tx"]]), ## Groups = E3 vs E4 cells
        name="Petropoulos et al's data",
        t2g=rowData(EMTAB3929[["tx"]])$gene
    ),
    list(
        xp=assays(GSE64016)$tx,
        grps=colData(GSE64016)$source_name_ch1=="single H1-Fucci cell sorted from G1 phase of the cell cycle only", ## Groups = G1 vs S phase of the cell cycle
        name="GSE64016",
        t2g=rowData(GSE64016[["tx"]])$gene
    )
)

## data = lapply(data,function(x){
##     x$xp = x$xp[1:5000,]
##     x$t2g = x$t2g[1:5000]
##     x
## })

## Choose threshold for inclusion of transcripts in the analysis (expressed in what fraction of the cells)
cutoffs_expressed=c(0.001,0.2)

## Compute TPMs, expression frequency of transcripts, binarized expression frequency, cellular detection rate
data=lapply(
    data,
    function(x){
        print(x$name)
        xp=x$xp
        dn=dimnames(xp)
        dn[[1]]=make.names(dn[[1]])
        if(is(xp,"sparseMatrix")){
            xp_log=xp
            xp_log@x=log2(1+xp_log@x)
            xp_tpm=xp%*%Diagonal(x=1/Matrix:::colSums(xp))*10^6 ## Compute TPMs
            xp_logtpm=xp_tpm
            xp_logtpm@x=log2(1+xp_logtpm@x)
        } else {
            xp_log=log2(1+xp)
            xp_tpm=xp%*%diag(1/colSums(xp))*10^6 ## Compute TPMs
            xp_logtpm=log2(1+xp_tpm)
        }
        expression_frequency=rowSums(as.matrix(xp)>0)/ncol(xp) ## Per transcript, in what fraction of the cells in the transcript expressed
        sca=SingleCellExperiment(
            assays=list(counts=as.matrix(xp),logcounts=as.matrix(xp_log),tpm=as.matrix(xp_tpm),logtpm=as.matrix(xp_logtpm)),
            colData=data.frame(
                Group=x$grps, ## Cell type labels
                CDR=colSums(as.matrix(xp)>0) ## Cellular Detection Rate (See Finak et al, Genome Biology, 2015)
            ),
            rowData=data.frame(
                rownames=dimnames(xp)[[1]],
                transcript=dimnames(xp)[[1]], ## Needed for do_mast_tx wrapper
                gene_ensembl=as.character(x$t2g), ## Transcript to gene mapping
                expression_frequency=expression_frequency,
                sapply(cutoffs_expressed,function(co){expression_frequency>co}), ## For each threshold, boolean indicating whether the transcript passes the threshold or not
                stringsAsFactors=FALSE                
            )
        )
        nc=ncol(rowData(sca))
        colnames(rowData(sca))[seq(nc-length(cutoffs_expressed)+1,nc,by=1)]=make.names(paste0("binary_",cutoffs_expressed))
        x=c(x,list(sca=sca))
    }
)
data=setNames(data,sapply(data,"[","name"))

rm(list=setdiff(ls(),c("data","do_logreg","do_mast_gn","do_mast_tx","do_t_tests","sidak","splatter_dir","cutoffs_expressed", "do_edgeRCDR", "do_edgeR", "do_limma", "do_DESeq2")))
gc()


####################
## 2. Perform GDE
####################

## For each dataset (loop over data)
##     For each threshold for transcript expression (loop over cutoffs_expressed)
##          Perform multivariate Log Reg
##          Perform MAST
##              For each GDE method, run it with or without CDR as a covariate, and with or without shuffling of the cell labels
##              Return vertors of p-values
file=file.path(splatter_dir,"shuffle_rds","degs.rds")
if(!file.exists(file)){
    data2=lapply(
        data,
        function(x){
            print(x$name)
            
            L=lapply(
                make.names(paste0("binary_",cutoffs_expressed)),
                function(co){
                    set.seed(123)
                    sca=x$sca
                    sca=sca[rowData(sca)[,co],]
                    ## Create a vector of shuffled cell types for future use
                    colData(sca)$Group_shuffled=sample(colData(sca)$Group)

                    ##
                    ## Below is logistic regression. It is similar to do_logreg() defined in wrappers.R
                    ##
                    xp_split=split_matrix(assays(sca)$counts,as.character(rowData(sca)$gene_ensembl)) ## split matrix into rows mapping to the same gene
                    xp_split=lapply(xp_split,t)
                    xp_split=lapply(xp_split,function(y){
                        colnames(y)=make.names(colnames(y)) ## Make sure that your column names are compatible with R formulas
                        y
                    })
                    xp_split=lapply(xp_split,cbind,grps=colData(sca)$Group,CDR=colData(sca)$CDR) ## Append CDR and cell types as new columns
                    xp_split=lapply(xp_split,as.data.frame)

                    lr_noCDR_noscrbl=unlist(mclapply(
                        xp_split,
                        function(y){
                            fmla=paste0("grps~",paste0(setdiff(colnames(y),c("grps","CDR")),collapse="+")) ## We remove CDR from the formula since we do not want to include it
                            lrtest(glm(formula=fmla,data=y,family=binomial()),glm("grps~1",data=y,family=binomial()))[2,5]
                        },
                        mc.cores=getOption("mc.cores",1L)
                    ))

                    ##
                    ## Below is MAST with or without CDR
                    ##
                    mast_CDR_noscrbl=do_mast_tx(
                        fmla="~Group+CDR",
                        sca=FromMatrix(assays(sca)$logtpm,cData=colData(sca),fData=rowData(sca)),
                        LRT="GroupTRUE"
                    )
                    
                    ## edgeRQLF_CDR
                    edgeRQLF_CDR = do_edgeRCDR(assays(sca)$counts, colData(sca)$Group)
                    edgeRQLF_CDR_noscrbl = list(
                        tx = edgeRQLF_CDR,
                        gn = sapply(
                            split(edgeRQLF_CDR, rowData(sca)[names(edgeRQLF_CDR),"gene_ensembl"]),
                            sidak
                        )
                    )
                    
                    ## limmatrend
                    limmatrend = do_limma(assays(sca)$counts, colData(sca)$Group)                
                    limmatrend_noscrbl = list(
                        tx = limmatrend,
                        gn = sapply(
                            split(limmatrend, rowData(sca)[names(limmatrend),"gene_ensembl"]),
                            sidak
                        )
                    )
                    
                    ##DESeq2
                    DESeq2 = do_DESeq2(assays(sca)$counts, colData(sca)$Group)                
                    DESeq2_noscrbl = list(
                        tx = DESeq2,
                        gn = sapply(
                            split(DESeq2, rowData(sca)[names(DESeq2),"gene_ensembl"]),
                            sidak
                        )
                    )

                    ##
                    ## Below is logistic regression, with shuffling of cell types.
                    ##
                    xp_split=split_matrix(assays(sca)$counts,as.character(rowData(sca)$gene_ensembl))
                    xp_split=lapply(xp_split,t)
                    xp_split=lapply(xp_split,function(y){
                        colnames(y)=make.names(colnames(y))
                        y
                    })
                    xp_split=lapply(xp_split,cbind,grps=colData(sca)$Group_shuffled,CDR=colData(sca)$CDR) ## We use Group_shuffled instead of Group
                    xp_split=lapply(xp_split,as.data.frame)

                    lr_noCDR_scrbl=unlist(mclapply(
                        xp_split,
                        function(y){
                            fmla=paste0("grps~",paste0(setdiff(colnames(y),c("grps","CDR")),collapse="+"))
                            lrtest(glm(formula=fmla,data=y,family=binomial()),glm("grps~1",data=y,family=binomial()))[2,5]
                        },
                        mc.cores=getOption("mc.cores",1L)
                    ))
                    
                    ##
                    ## Below is MAST, with shuffling of cell types.
                    ##
                    mast_CDR_scrbl=do_mast_tx(
                        fmla="~Group_shuffled+CDR",
                        sca=FromMatrix(assays(sca)$logtpm,cData=colData(sca),fData=rowData(sca)),
                        LRT="Group_shuffledTRUE"
                    )

                    ## edgeRQLF_CDR
                    edgeRQLF_CDR = do_edgeRCDR(assays(sca)$counts, colData(sca)$Group_shuffled)
                    edgeRQLF_CDR_scrbl = list(
                        tx = edgeRQLF_CDR,
                        gn = sapply(
                            split(edgeRQLF_CDR, rowData(sca)[names(edgeRQLF_CDR),"gene_ensembl"]),
                            sidak
                        )
                    )
                    
                    ## limmatrend
                    limmatrend = do_limma(assays(sca)$counts, colData(sca)$Group_shuffled)                
                    limmatrend_scrbl = list(
                        tx = limmatrend,
                        gn = sapply(
                            split(limmatrend, rowData(sca)[names(limmatrend),"gene_ensembl"]),
                            sidak
                        )
                    )
                    
                    ##DESeq2
                    DESeq2 = do_DESeq2(assays(sca)$counts, colData(sca)$Group_shuffled)
                    DESeq2_scrbl = list(
                        tx = DESeq2,
                        gn = sapply(
                            split(DESeq2, rowData(sca)[names(DESeq2),"gene_ensembl"]),
                            sidak
                        )
                    )

                    return_names=apply(expand.grid(c("mast_CDR","lr_noCDR", "edgeRQLF_CDR", "limmatrend", "DESeq2"),c("noscrbl","scrbl")),1,paste0,collapse="_")
                    
                    return(mget(return_names))
                }
            )
            names(L)=cutoffs_expressed
            x=c(x,degs=list(L))
            x
        }
    )
    saveRDS(data2,file=file)
    data=data2
} else {
    data=readRDS(file)
}

####################
## 2. Perform GDE on cells from the same batch (randomly affected to two distinct subgroups)
####################

## For each dataset (loop over data)
##     For each threshold for transcript expression (loop over cutoffs_expressed)
##         For each of the two cell population
##              Subset on the cell population
##              Create a random vector of cell type
##              Perform multivariate Log Reg (trying to predict a random boolean vector)
##              Perform MAST  (trying to predict a random boolean vector)
##                  For each GDE method, run it with or without CDR as a covariate
##                  Return vertors of p-values
file=file.path(splatter_dir,"shuffle_rds","neg.rds")
if(!file.exists(file)){
    data2=lapply(
        data,
        function(x){
            print(x$name)
            L=lapply(
                make.names(paste0("binary_",cutoffs_expressed)),
                function(co){
                    set.seed(123)
                    sca=x$sca
                    sca=sca[rowData(sca)[,co],]
                    sapply(
                        unique(colData(sca)$Group),
                        function(celltype){
                            ## Subset on cell types
                            sca=sca[,colData(sca)$Group==celltype]

                            ## Create vector with 50% TRUE and 50% FALSE, randomly ordered
                            w=rep(FALSE,ncol(sca))
                            n=ceiling(length(w)/2)
                            w[sample(1:length(w),n)]=TRUE
                            response=as.numeric(w)
                            colData(sca)$Group_shuffled=w

                            ##
                            ## Run logistic regression, same as above
                            ##
                            xp_split=split_matrix(assays(sca)$counts,as.character(rowData(sca)$gene_ensembl))
                            xp_split=lapply(xp_split,t)
                            xp_split=lapply(xp_split,function(y){
                                colnames(y)=make.names(colnames(y))
                                y
                            })
                            xp_split=lapply(xp_split,cbind,grps=colData(sca)$Group_shuffled,CDR=colData(sca)$CDR)
                            xp_split=lapply(xp_split,as.data.frame)

                            lr_noCDR_scrbl=unlist(mclapply(
                                xp_split,
                                function(y){
                                    fmla=paste0("grps~",paste0(setdiff(colnames(y),c("grps","CDR")),collapse="+"))
                                    lrtest(glm(formula=fmla,data=y,family=binomial()),glm("grps~1",data=y,family=binomial()))[2,5]
                                },
                                mc.cores=getOption("mc.cores",1L)
                            ))

                            
                            ##
                            ## Run MAST, same as above
                            ##
                            mast_CDR_scrbl=do_mast_tx(
                                fmla="~Group_shuffled+CDR",
                                sca=FromMatrix(assays(sca)$logtpm,cData=colData(sca),fData=rowData(sca)),
                                LRT="Group_shuffledTRUE"
                            )
                            
                            ## edgeRQLF_CDR
                            edgeRQLF_CDR = do_edgeRCDR(assays(sca)$counts, colData(sca)$Group_shuffled)
                            edgeRQLF_CDR_scrbl = list(
                                tx = edgeRQLF_CDR,
                                gn = sapply(
                                    split(edgeRQLF_CDR, rowData(sca)[names(edgeRQLF_CDR),"gene_ensembl"]),
                                    sidak
                                )
                            )
                            
                            ## limmatrend
                            limmatrend = do_limma(assays(sca)$counts, colData(sca)$Group_shuffled)                
                            limmatrend_scrbl = list(
                                tx = limmatrend,
                                gn = sapply(
                                    split(limmatrend, rowData(sca)[names(limmatrend),"gene_ensembl"]),
                                    sidak
                                )
                            )
                            
                            ##DESeq2
                            DESeq2 = do_DESeq2(assays(sca)$counts, colData(sca)$Group_shuffled)                
                            DESeq2_scrbl = list(
                                tx = DESeq2,
                                gn = sapply(
                                    split(DESeq2, rowData(sca)[names(DESeq2),"gene_ensembl"]),
                                    sidak
                                )
                            )
                                                       
                            return_names=apply(expand.grid(c("mast_CDR","lr_noCDR", "edgeRQLF_CDR", "limmatrend", "DESeq2"),c("scrbl")),1,paste0,collapse="_")
                            
                            ##return(mget(return_names))
                            return(list(pvalues=mget(return_names),scrbl=colData(sca)[,"Group_shuffled",drop=FALSE]))
                        },
                        simplify=FALSE
                    )
                }
            )
            names(L)=cutoffs_expressed
            x=c(x,negs=list(L))
            x
        }
    )
    saveRDS(data2,file=file)
    data=data2
} else {
    data=readRDS(file)
}
## Appending names of celltypes
data=lapply(
    data,
    function(x){
        groups=unique(colData(x$sca)$Group)
        x$negs=lapply(x$negs,setNames,groups)
        x
    }
)


####################
## 4. Plots
####################

if(FALSE){
    ##
    ## Power analysis (Freq. of discovery vs nominal FDR)
    ##
    
    methods=apply(expand.grid(c("mast_CDR","lr", "edgeR", "limma", "DESeq2"),c("noscrbl","scrbl")),1,paste0,collapse="_")
    colors=setNames(rep(brewer.pal(5, "Set1"),2),methods)
    ## linetype=setNames(rep(c("solid","dotted","dashed","dotdash")[1:2],each=length(methods)/2),methods)
    ## linetype = setNames(ifelse(grepl("_noscrbl", methods,), "solid", "dotted"), methods)
    linetype = setNames(ifelse(grepl("_noscrbl", methods,), "solid", "dotted"), methods)
    ## Loop over GDE results, extract p value vectors
    ##     Compute freq of discovery for each p-value threshold
    ##     Return a ggplot object(ggplot2) for each dataset, expression cutoff
    for(cutoff in as.character(cutoffs_expressed)){
        file=file.path(splatter_dir,"graphs",paste0("power_analysis_co",cutoff,".png"))
        png(filename=file,height=3000,width=3000,res=300)
        ps=lapply(names(data),function(ds){
            ## Extract p value vectors
            df=data.table(
                lr_noscrbl=data[[ds]]$degs[[cutoff]]$lr_noCDR_noscrbl,
                ## lr_CDR_noscrbl=data[[ds]]$degs[[cutoff]]$lr_CDR_noscrbl$gn,
                lr_scrbl=data[[ds]]$degs[[cutoff]]$lr_noCDR_scrbl,
                ## lr_CDR_scrbl=data[[ds]]$degs[[cutoff]]$lr_CDR_scrbl$gn,
                ## mast_noCDR_noscrbl=data[[ds]]$degs[[cutoff]]$mast_noCDR_noscrbl$gn,
                mast_CDR_noscrbl=data[[ds]]$degs[[cutoff]]$mast_CDR_noscrbl$gn,
                ## mast_noCDR_scrbl=data[[ds]]$degs[[cutoff]]$mast_noCDR_scrbl$gn,
                mast_CDR_scrbl=data[[ds]]$degs[[cutoff]]$mast_CDR_scrbl$gn,
                edgeR_noscrbl=data[[ds]]$degs[[cutoff]]$edgeRQLF_CDR_noscrbl$gn,
                edgeR_scrbl=data[[ds]]$degs[[cutoff]]$edgeRQLF_CDR_scrbl$gn,
                limma_noscrbl=data[[ds]]$degs[[cutoff]]$limmatrend_noscrbl$gn,
                limma_scrbl=data[[ds]]$degs[[cutoff]]$limmatrend_scrbl$gn,
                DESeq2_noscrbl=data[[ds]]$degs[[cutoff]]$DESeq2_noscrbl$gn,
                DESeq2_scrbl=data[[ds]]$degs[[cutoff]]$DESeq2_scrbl$gn
            )
            ## Compute freq of discovery
            df=suppressWarnings(melt.data.table(df))
            colnames(df)=c("method","p")
            df=split(df,df$method)
            df=lapply(df,function(x){
                x[is.nan(p),p:=1]
                x[,o:=order(p)]
                x[,q:=p.adjust(p,method="fdr")]
                x[,n:=nrow(.SD),by=p]
                setorder(x,p)
                x=x[!duplicated(p),]
                x=x[,n:=cumsum(x$n)]
                x=x[,n:=n/max(n)]
                x
            })    
            df=do.call(rbind,df)

            ## Generate plot objects
            p1=ggplot(data=df,aes(x=p, y=n, colour=method, linetype=method)) +
                geom_abline(slope=max(df$n),intercept=0) +
                geom_path(size=1) +
                scale_colour_manual(values=colors) +
                scale_linetype_manual(values=linetype) +
                xlab("p") +
                ylab("# discoveries")+
                theme_bw() +
                ggtitle(ds)

            p2=ggplot(data=df,aes(x=q, y=n, colour=method, linetype=method)) +
                geom_abline(slope=max(df$n),intercept=0) +
                geom_vline(xintercept=0.05,linetype="dotted") +
                geom_path(size=1) +
                scale_colour_manual(values=colors) +
                scale_linetype_manual(values=linetype) +
                xlab("q") +
                ylab("# discoveries")+
                theme_bw() +
                ggtitle("")
            list(p1,p2)
        })
        ps=do.call(c,ps)
        print(ggarrange(plotlist=ps,ncol=2,nrow=3))
        dev.off()
    }

    ## Power analysis for shuffled data. Similar to above expect we use a different set of pvalue vectors
    for(cutoff in as.character(cutoffs_expressed)){
        file=file.path(splatter_dir,"graphs",paste0("neg_power_analysis_co",cutoff,".png"))
        png(filename=file,height=3000,width=6000,res=300)
        ps=lapply(names(data),function(ds){
            do.call(c,lapply(as.character(c(TRUE,FALSE)),function(pop){
                ## Only difference from above is the source of data parsing (and number of methods)
                df=data.table(
                    lr_noCDR_scrbl=data[[ds]]$negs[[cutoff]][[pop]]$pvalues$lr_noCDR_scrbl$gn,
                    lr_CDR_scrbl=data[[ds]]$negs[[cutoff]][[pop]]$pvalues$lr_CDR_scrbl$gn,
                    mast_noCDR_scrbl=data[[ds]]$negs[[cutoff]][[pop]]$pvalues$mast_noCDR_scrbl$gn,
                    mast_CDR_scrbl=data[[ds]]$negs[[cutoff]][[pop]]$pvalues$mast_CDR_scrbl$gn        
                )
                df=suppressWarnings(melt.data.table(df))
                colnames(df)=c("method","p")
                df=split(df,df$method)
                df=lapply(df,function(x){
                    x[is.nan(p),p:=1]
                    x[,o:=order(p)]
                    x[,q:=p.adjust(p,method="fdr")]
                    x[,n:=nrow(.SD),by=p]
                    setorder(x,p)
                    x=x[!duplicated(p),]
                    x=x[,n:=cumsum(x$n)]
                    x=x[,n:=n/max(n)]
                    x
                })    
                df=do.call(rbind,df)
                
                p1=ggplot(data=df,aes(x=p, y=n, colour=method, linetype=method)) +
                    geom_abline(slope=max(df$n),intercept=0) +
                    geom_path(size=1) +
                    scale_colour_manual(values=colors) +
                    scale_linetype_manual(values=linetype) +
                    xlab("p") +
                    ylab("# discoveries")+
                    theme_bw() +
                    annotate("text", -Inf, Inf, label = ifelse(pop,"CellPop1","CellPop2"), hjust = 0, vjust = 1)
                ggtitle(ifelse(pop,ds,""))

                p2=ggplot(data=df,aes(x=q, y=n, colour=method, linetype=method)) +
                    geom_abline(slope=max(df$n),intercept=0) +
                    geom_vline(xintercept=0.05,linetype="dotted") +
                    geom_path(size=1) +
                    scale_colour_manual(values=colors) +
                    scale_linetype_manual(values=linetype) +
                    xlab("q") +
                    ylab("# discoveries")+
                    theme_bw() +
                    ggtitle("")
                list(p1,p2)
            }))
        })
        ps=do.call(c,ps)
        print(ggarrange(plotlist=ps,ncol=4,nrow=3))
        dev.off()
    }

    ##
    ## CDR histograms
    ##

    ## Plot distribution of CDRs as histograms for each dataset
    file=file.path(splatter_dir,"graphs",paste0("CDR_histograms.png"))
    png(file,width=3000,height=1000,res=300)
    par(mfrow=c(1,3))
    sapply(data,function(x){
        d=colData(x$sca)$CDR/nrow(x$sca)
        hist(d,main=x$name,xlab="CDR",col="lightgray",freq=FALSE,xlim=c(0,max(d)))
    })
    dev.off()

    ##
    ## Power analysis stratified by expression frequency
    ##

    ## Here we study each method's output on univariate transcripts (or TCCs), binned by expression frequency
    cutoff="0.001"
    breaks=c(0.001,0.01,0.1,0.2,1) ## Define bins boundaries
    n_strata=length(breaks)-1

    ## Define color, linetypes...
    colors=colorRampPalette(c("dodgerblue","chartreuse3","darkorange","brown"))(n_strata)
    linetype=setNames(rep(c("solid","dotted","dashed","dotdash"),each=4),methods)

    file=file.path(splatter_dir,"graphs",paste0("power_analysis_stratified.png"))
    png(filename=file,height=3500,width=7000,res=300)
    ps=lapply(names(data),function(ds){
        print(ds)
        ## Load DEG results
        df=data.table(
            lr_noCDR_noscrbl=data[[ds]]$degs[[cutoff]]$lr_noCDR_noscrbl$tx,
            ## lr_CDR_noscrbl=data[[ds]]$degs[[cutoff]]$lr_CDR_noscrbl$tx,
            ## mast_noCDR_noscrbl=data[[ds]]$degs[[cutoff]]$mast_noCDR_noscrbl$tx,
            mast_CDR_noscrbl=data[[ds]]$degs[[cutoff]]$mast_CDR_noscrbl$tx
        )

        ## Load expression frequencies and bin it
        freqs=rowData(data[[ds]]$sca)$expression_frequency[rowData(data[[ds]]$sca)[,paste0("binary_",cutoff)]]
        freqs=cut(freqs,breaks=breaks,include.lowest=TRUE)
        names(colors)=levels(freqs)
        df=cbind(df,freqs)

        ## Compute frequency of discovery for each possible p-value threshold
        df=suppressWarnings(melt.data.table(df))
        colnames(df)=c("freqs_bin","method","p")
        df$group=paste0(df$method,df$freqs_bin)
        df$alg=sapply(strsplit(df$group,"_"),"[",1)
        df=split(df,df$group)
        df=lapply(df,function(x){
            x[is.nan(p),p:=1]
            x[,o:=order(p)]
            x[,q:=p.adjust(p,method="fdr")]
            x[,n:=nrow(.SD),by=p]
            setorder(x,p)
            x=x[!duplicated(p),]
            x=x[,n:=cumsum(x$n)]
            x[,n:=n/max(n)]
            x
        })    
        df=do.call(rbind,df)

        ## Generate power analysis plots
        p1=ggplot(data=subset(df,grepl("mast",df$alg)),aes(x=p, y=n, group=group,colour=freqs_bin, linetype=method)) +
            geom_path(size=1) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            xlab("p") +
            ylab("# discoveries")+
            theme_bw() +
            ggtitle(ds) +
            xlim(0,1) +
            annotate("text", -Inf, Inf, label = "MAST", hjust = 0, vjust = 1)

        p2=ggplot(data=subset(df,grepl("lr",df$alg)),aes(x=p, y=n, group=group,colour=freqs_bin, linetype=method)) +
            geom_path(size=1) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            xlab("p") +
            ylab("# discoveries")+
            theme_bw() +
            xlim(0,1) +
            annotate("text", -Inf, Inf, label = "LR", hjust = 0, vjust = 1)

        
        p3=ggplot(data=subset(df,grepl("mast",df$alg)),aes(x=q, y=n, group=group, colour=freqs_bin, linetype=method)) +
            geom_vline(xintercept=0.05,linetype="dotted") +
            geom_path(size=1) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            xlab("q") +
            ylab("# discoveries")+
            theme_bw() +
            ggtitle("") +
            xlim(0,1) +
            annotate("text", -Inf, Inf, label = "MAST", hjust = 0, vjust = 1)


        p4=ggplot(data=subset(df,grepl("lr",df$alg)),aes(x=q, y=n, group=group, colour=freqs_bin, linetype=method)) +
            geom_vline(xintercept=0.05,linetype="dotted") +
            geom_path(size=1) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            xlab("q") +
            ylab("# discoveries")+
            theme_bw() +
            ggtitle("") +
            xlim(0,1) +
            annotate("text", -Inf, Inf, label = "LR", hjust = 0, vjust = 1)

        list(p1,p2,p3,p4)
    })
    ps=do.call(c,ps)
    print(ggarrange(plotlist=ps,ncol=4,nrow=3))
    dev.off()

    ##
    ## Histograms of p-values for the random cell labels
    ##
    DT=data.table()
    for(ds in names(data)){
        print(ds)
        for(co in as.character(cutoffs_expressed)){
            for(celltype in names(data[[ds]]$negs[[co]])){
                for(feature_type in c("tx","gn")){
                    dt=as.data.table(sapply(data[[ds]]$negs[[co]][[celltype]]$pvalues,function(x){x[[feature_type]]}))
                    dt=suppressWarnings(melt(dt))
                    colnames(dt)=c("method","p")
                    dt[,c("celltype","co","ds"):=list(celltype,co,ds)]
                    DT=rbind(DT,dt)
                }
            }
        }
    }

    file=file.path(splatter_dir,"graphs",paste0("neg_histograms.pdf"))
    pdf(file,width=14)
    invisible(apply(
        unique(DT[,c("celltype","ds")]),
        1,
        function(slice){
            p=ggplot(data=DT[celltype==slice["celltype"]&ds==slice["ds"],],aes(x=p)) +
                geom_histogram() +
                facet_wrap(~ds+co+method,scale="free_y",ncol=length(unique(DT$method)),nrow=length(unique(DT$fold))) +
                theme_bw() +
                xlab("p") +
                ylab("Counts per bin") +
                xlim(0,1) +
                ggtitle(ifelse(slice["celltype"]=="TRUE","CellPop1","CellPop2"))
            print(p)
        }
    ))
    dev.off()
}
