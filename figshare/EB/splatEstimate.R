{
    library(splatter)
    library(MultiAssayExperiment)
    library(rio)
    library(biomaRt)
    library(ebmisc)
    library(Matrix)
    library(MAST)
    library(ggplot2)
    library(scater)
    library(ebmisc)
    library(splatter)
    library(gridExtra)
    library(reshape)
    library(grid)
    library(RColorBrewer)
    library(ggpubr)
    library(patchwork)
    library(tools)
    library(TENxPBMCData)

    options(mc.cores=8L)

    ## Load wrappers to run LR & MAST
    source("../../wrappers.R")

    ## Adapted from: https://github.com/pachterlab/NYMP_2018/blob/master/simulations/RSEM/R/roc_helpers.R
    source("../../NYMP_2018/simulations/RSEM/R/roc_helpers.R")

    ## Get transcript to gene mapping
    splatter_dir="."
    t2g=readRDS(file=file.path("../splatter/","t2g.rds"))
    rownames(t2g)=t2g$ensembl_transcript_id

    params_dir="./rds"
    dir.create(params_dir,showWarnings=FALSE)
    dir.create(file.path(splatter_dir,"input"),showWarnings=FALSE)
    dir.create(file.path(splatter_dir,"pre_input"),showWarnings=FALSE)

    if(FALSE){
        library(TENxPBMCData)
        library(TabulaMurisData)
        suppressPackageStartupMessages({
            library(ExperimentHub)
            library(SingleCellExperiment)
            library(TabulaMurisData)
        })
        
        TENxPBMC4k=TENxPBMCData("pbmc4k")
        saveRDS(as.matrix(counts(TENxPBMC4k)),file=file.path(splatter_dir,"input","10xv2PBMC4k.Rds"))

        droplet <- TabulaMurisDroplet()
        organ="Thymus"
        data <- droplet[, droplet$tissue == "Thymus"]
        counts <- as.matrix(counts(data))
        saveRDS(counts,file=file.path(splatter_dir,"input","TM_thymus.Rds"))

	## Requires Seurat v3
        library(Seurat)
        data1k_v3 <- Read10X(file.path(splatter_dir,"pre_input/pbmc_1k_v3_filtered_feature_bc_matrix/filtered_feature_bc_matrix")) ## Source: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3
        counts <- as.matrix(data1k_v3)
        saveRDS(counts,file=file.path(splatter_dir,"input","10xv3PBMC1k.Rds"))
    }
    
    files=list.files("./input",full.names=TRUE)
    if(FALSE){
        for(file in files){

            print(file)

            d=readRDS(file)

            set.seed(123)
            params=newSplatParams()
            params=try(splatEstimate(d,params))

            saveRDS(
                params,
                file=file.path(params_dir,paste0("splatEstimate_",basename(file)))
            )
            
        }
    }

    nGenes=nrow(t2g)
    params=lapply(
        files,
        function(file){
            param=readRDS(file.path(params_dir,paste0("splatEstimate_",basename(file))))
            param=setParams(
                param,
                nGenes=nGenes,
                batchCells=500,
                group.prob=c(0.5,0.5),
                de.facLoc=1.5
            )
            list(name=sub(".Rds","",basename(file),fixed=TRUE),params=param)
        } 
    )
    t2g=t2g[1:nGenes,]

    simulations=lapply(params,function(x){
        splat=splatSimulate(x$params,seed=123,method="group")
        rownames(splat)=t2g$ensembl_transcript_id
        x=c(x,list(splat=splat))
        x
    })
    names(simulations)=sapply(simulations,function(x){x$name})

    ##
    ## Computing tpms, cntxon...
    ##
    simulations=lapply(
        simulations,function(x){
            splat=x$splat
            colData(splat)$cntxon=colSums(assays(splat)$counts>0) ##cntxon = Cellular Detection Rate (CDR) at the transcript level (see Finak et al, Genome Biology, 2015)
            rowData(splat)$transcript=rownames(rowData(splat)) ## Necessary for the do_mast_tx wrapper
            rowData(splat)$gene_ensembl=t2g[rownames(rowData(splat)),"ensembl_gene_id"] ## Necessary for the do_mast_tx wrapper
            assays(splat)$tpm=assays(splat)$count%*%diag(1/colSums(assays(splat)$count))*10^6 ## Computing TPMs
            rowData(splat)$expr_freq=rowMeans(assays(splat)$count>0)
            splat=splat[rowData(splat)$expr_freq>0,] ## Remove genes with only 0s
            x$splat=splat
            x
        }
    )

####################
    ## 4. Add gene-aggregated data
####################

    simulations=lapply(
        simulations,
        function(x){
            count=aggregate_df(assays(x$splat)$count,rowData(x$splat)$gene_ensembl,fun=sum,margin=1) ## Aggregate transcripts to gene counts by summation
            cData=colData(x$splat)
            fData=rowData(x$splat)
            fData=aggregate_df(as.matrix(fData)[,grepl("DEFacGroup",colnames(fData))],fData$gene_ensembl,fun=function(y){
                res=apply(!as.matrix(y)==1,2,any)
                mode(res)="numeric"
                res
            }) ## rowData only keeps track of whether the gene is differentially expressed (it is if any of its transcript is)
            fData=cbind(fData,gene_ensembl=rownames(fData),transcript=rownames(fData)) ## Necessary for the do_mast_tx wrapper
            splat_gn=FromMatrix(count,cData=cData,fData=as(fData,"DataFrame"),check_sanity=FALSE) ## Create the SingleCellAssay object
            names(assays(splat_gn))="counts"
            colData(splat_gn)$cntxon=colSums(assays(splat_gn)$counts>0) ##cntxon = Cellular Detection Rate (CDR) at the transcript level (see Finak et al, Genome Biology, 2015)
            rowData(splat_gn)$transcript=rownames(rowData(splat_gn))
            rowData(splat_gn)$gene_ensembl=rownames(rowData(splat_gn))
            assays(splat_gn)$tpm=assays(splat_gn)$count%*%diag(1/colSums(assays(splat_gn)$count))*10^6 ## compute TPMs
            x$splat_gn=splat_gn
            x
        }
    )

    ## Export all simulations
    ## Export parameters
    dir.create(file.path(splatter_dir,"tables"), showWarnings=FALSE)
    write.table(
        do.call(
            rbind,
            sapply(
                simulations,
                function(x){
                    data.frame(
                        name=x$name,
                        Group1Batch1=nrow(subset(colData(x$splat),Group=="Group1"&Batch=="Batch1")),
                        Group2Batch1=nrow(subset(colData(x$splat),Group=="Group2"&Batch=="Batch1")),
                        dropout.type=x$params@dropout.type,
                        dropout.midrange=paste0(signif(x$params@dropout.mid,2),collapse=" / "),
                        dropout.shape=paste0(signif(x$params@dropout.shape,2),collapse=" / "),
                        mean.shape=x$params@mean.shape,
                        mean.rate=x$params@mean.rate,
                        DE.facLoc=paste0(x$params@de.facLoc,collapse=" / "),
                        bcv.common=x$params@bcv.common,
                        bcv.df=x$params@bcv.df,
                        out.prob=x$params@out.prob,
                        out.facLoc=x$params@out.facLoc,
                        out.facScale=x$params@out.facScale
                    )
                },
                simplify=FALSE
            )
        ),
        file=file.path(splatter_dir,"tables","simulations_parameters.csv"),
        sep=",",
        row.names=FALSE
    )

####################
    ## 5. Perform DEGs
####################
    simulations=lapply(
        simulations,
        function(x){
            print(x$name)
            ## Only perform DEGs if result is missing from disk, otherwise load.
            file=file.path(splatter_dir,"rds",paste0(x$name,"_degs.rds"))
            if(!file.exists(file)){
            ##if(TRUE){
                ##
                ## MAST
                ##
                fmla=if(length(unique(colData(x$splat)$Batch))>1){"~Group+Batch+cntxon"}else{"~Group+cntxon"} ## Model is (transcript counts ~ Group + Batch) if we simulate Batch, otherwise it is just (transcript counts ~ Group). We do not use cntxon as it is not informative in splatter so there is no reason to include it in the model
             
                sca=FromMatrix(log2(1+assays(x$splat)$tpm),cData=colData(x$splat),fData=rowData(x$splat))
                LRT="GroupGroup2"
                zlmCond=zlm(as.formula(fmla),sca)
                mast_tpm=lrTest(zlmCond,CoefficientHypothesis("GroupGroup2"))
                mast_tpm=mast_tpm[,"hurdle","Pr(>Chisq)"]
                mast_tpm=list(
                    gn=unlist(sapply(
                        split(mast_tpm,rowData(x$splat)[names(mast_tpm),"gene_ensembl"]),
                        function(slice){
                            ## print(slice)
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
                    ),use.names=FALSE),
                    tx=mast_tpm
                )
                mast_tpm=lapply(mast_tpm,function(x){x[is.na(x)]=1;x})
                
                sca=FromMatrix(log2(1+assays(x$splat_gn)$tpm),cData=colData(x$splat_gn),fData=rowData(x$splat_gn))
                LRT="GroupGroup2"
                zlmCond=zlm(as.formula(fmla),sca)
                mast_unaggregated_gn_tpm=lrTest(zlmCond,CoefficientHypothesis("GroupGroup2"))
                mast_unaggregated_gn_tpm=mast_unaggregated_gn_tpm[,"hurdle","Pr(>Chisq)"]
                mast_unaggregated_gn_tpm[is.na(mast_unaggregated_gn_tpm)]=1
                mast_unaggregated_gn_tpm=list(gn=mast_unaggregated_gn_tpm)
                
                rm(zlmCond,sca,LRT)
                

                ## Use Batch as covariate if relevant, NULL otherwise
                if(length(unique(colData(x$splat)[,"Batch",drop=TRUE]))>1){
                    covariates=colData(x$splat)[,c("Batch","cntxon"),drop=FALSE]
                    covariates_aov=covariates[,1,drop=TRUE]
                } else {
                    covariates=colData(x$splat)[,c("cntxon"),drop=FALSE]
                    covariates_aov=rep(1,ncol(x$splat))
                }

                
                ##
                ## multivariate Logistic Regression
                ##
                lr_tpm=do_logreg(data=1+assays(x$splat)$tpm,response=colData(x$splat)$Group=="Group2",tx2gene=rowData(x$splat)$gene_ensembl,covariates=covariates) ## mLR on transcripts
                lr_unaggregated_tpm=do_logreg(data=1+assays(x$splat)$tpm,response=colData(x$splat)$Group=="Group2",tx2gene=rowData(x$splat)$transcript,covariates=covariates) ## univariate LR on transcripts
                lr_unaggregated_gn_tpm=do_logreg(data=1+assays(x$splat_gn)$tpm,response=colData(x$splat_gn)$Group=="Group2",tx2gene=rowData(x$splat_gn)$transcript,covariates=covariates) ## unvariate LR on gene counts

                ## ANOVA
                aov_tpm=do_t_tests(response=colData(x$splat)$Group,data=log2(1+assays(x$splat)$tpm),t2g=rowData(x$splat)$gene_ensembl,covariates=covariates_aov)
                aov_unaggregated_gn_tpm=do_t_tests(response=colData(x$splat_gn)$Group,data=log2(1+assays(x$splat_gn)$tpm),t2g=rowData(x$splat_gn)$gene_ensembl,covariates=covariates_aov)
                
                rm(list=c("fmla","covariates","covariates_aov"))
                env=environment()
                ls=setdiff(ls(env=env),c("env","file","x"))
                degs=sapply(ls,get,env=env,simplify=FALSE)
                x=c(x,list(degs=degs))
                saveRDS(x,file=file)
                x
            } else {
                readRDS(file)
            }
        }
    )

    ##
    ## Reformat DEGs (for consistency)
    ## 
    simulations=lapply(
        simulations,
        function(x){
            degs=x$degs
            x$degs_formatted=list(
                mast_tx_sidak=degs$mast_tpm$gn,
                mast_tx_univ=degs$mast_tpm$tx,
                mast_gn_univ=degs$mast_unaggregated_gn_tpm$gn,
                lr_tx_multiv=degs$lr_tpm,
                lr_tx_univ=degs$lr_unaggregated_tpm,
                lr_gn_univ=degs$lr_unaggregated_gn_tpm,
                aov_tx_sidak=degs$aov_tpm$gn,
                aov_tx_univ=degs$aov_tpm$tx,
                aov_gn_univ=degs$aov_unaggregated_gn_tpm$gn
            )
            x
        }
    )

    saveRDS(simulations,file=file.path(splatter_dir,"rds","simulations.rds"))

####################
    ## 6. Plots
####################

    lapply(
        simulations,
        function(x){
            
            sapply(file.path(splatter_dir,"graphs_paper",c("","genes_violin_plots","pvalues_histogram","roc_curves")),dir.create,showWarnings=FALSE)
            
            ## Compute which genes are differentially expressed (legacy, now also available in rowData(x$splat))
            print(x$name)
            de_transcripts=rowSums(as.matrix(rowData(x$splat)[,grepl("DEFacGroup",colnames(rowData(x$splat)))])!=1)>0
            de_genes=sapply(
                split_matrix(as.matrix(rowData(x$splat)),rowData(x$splat)$gene_ensembl),
                function(y){
                    any(rowSums(as.matrix(y[,grepl("DEFacGroup",colnames(y))])!=1)>0)
                }
            )
            rowData(x$splat)=cbind(rowData(x$splat),de_transcript=de_transcripts)
            rowData(x$splat)=cbind(rowData(x$splat),de_gene=de_genes[rowData(x$splat)$gene_ensembl])
            rowData(x$splat_gn)=cbind(rowData(x$splat_gn),de_genes=de_genes[rowData(x$splat_gn)$gene_ensembl])

            ## Generate 5-panel figures with 1) empirical FDR vs TPR [as per Nnatros et al], 2) FPR vs TPR [ROC curve], 3) nominal FDR (p-values from each GDE method + p.adjust(.,method="fdr") versus empirical FDR, 4) nominal FDR vs FPR and 5) nominal FDR vs sensitivity
            file=file.path(splatter_dir,"graphs_paper","roc_curves",paste0(x$name,".svg"))
            library(RSvgDevice)
            devSVG(file,height=3.5,width=3.5*2)
            fdrs=sapply(c("mast_tx_univ","mast_tx_sidak","lr_tx_univ","lr_tx_multiv"),function(y){
                if(grepl("ENSG",names(x$degs_formatted[[y]])[1])){
                    true_DE=de_genes
                } else if (grepl("ENST",names(x$degs_formatted[[y]])[1])){
                    true_DE=de_transcripts
                }
                calculate_fdr(true_DE=true_DE[names(x$degs_formatted[[y]])],pvalues=x$degs_formatted[[y]],title=y)
            },simplify=FALSE)

            fdrs <- lapply(1:length(fdrs), function(y) data.frame(sensitivity=fdrs[[y]]$sen, actual_FDR=fdrs[[y]]$fdr, FalsePositiveRate = fdrs[[y]]$FalsePositiveRate, nominal_FDR=p.adjust(fdrs[[y]]$pvalues,method="fdr"),Method = names(fdrs)[y]))
            fdrs <- do.call(rbind, fdrs)
            colors=setNames(makeTransparent(rep(brewer.pal(3,"Set1"),each=2),123),unique(fdrs$Method))
            linetype=setNames(rep(c("solid","dashed","dotted")[1:2],3),unique(fdrs$Method))
            p2 <- ggplot(data = fdrs, aes(x = FalsePositiveRate , y = sensitivity, colour=Method, linetype=Method)) +
                scale_colour_manual(values=colors) +
                scale_linetype_manual(values=linetype) +
                scale_x_sqrt() +
                geom_path() +
                geom_abline(intercept=0,slope=1,linetype="solid",col="gray") +
                theme_bw(12) +
                coord_fixed() +
                ggtitle(x$name)
            print(p2)
            dev.off()
        }
    )
}
