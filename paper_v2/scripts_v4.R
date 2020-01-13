library(patchwork)
library(gplots)
library(data.table)
library(ebmisc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gplots)
library(tools)
library(grid)
source("../wrappers.R")

options("mc.cores"=8L)

## Add transparency to a color
tp=function(hexcol,tp="44"){
    setNames(paste0(substr(hexcol,1,7),tp),names(hexcol))
}

####################
## Panels A), B) and C) CD45 differential analyses
####################

## Load results of the CD45 analysis, based on a script adapted from https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
load("./NYMP_2018/10x_example-logR/rdata/Fig2.RData")
load("./NYMP_2018/10x_example-logR/rdata/Fig2_allcells.RData")

## Reshape data
ps=lapply(ps,t)
ps=lapply(ps,as.data.table)
ps=suppressWarnings(lapply(ps,melt))
ps=sapply(names(ps),function(x){
    colnames(ps[[x]])=c("method","p")
    ps[[x]]$method=as.character(ps[[x]]$method)
    ps[[x]][,"celltypes":=x]
    ps[[x]][!grepl("_tx",method),]
},simplify=FALSE)

## Setup colors
cols=sapply(c("logreg_gene"="red","logreg_indep"="blue","logreg_tcc"="orange","mast_tcc"="purple"),col2hex)
cols=setNames(tp(cols,"99"),names(cols))

pdf("./graphs/PanelCD45.pdf",width=7,height=7/3)
par(mfrow=c(1,3),mar=c(3,3,3,1),cex.axis=0.9,cex.lab=1.2)
lapply(
    ## For each pairwise comparison of cell types
    names(ps),
    function(i){
        fig2a=copy(ps[[i]]) ## data.table uses assignment by reference, so we copy that part of the data first.
        fig2a$p=-log10(fig2a$p)
        xlim=c(0,ceiling(max(fig2a$p)))
        breaks=seq(xlim[1],xlim[2],length.out=20) ## For histograms, use common breaks
        
        fig2a=split(fig2a$p,fig2a$method)

        ## For each GDE method, compute the density of -log10(pvalues)
        densities=lapply(fig2a,density,from=0,to=ceiling(xlim[2]),n=1024)

        ## For each GDE method, plot histogram of pvalues
        hists=sapply(
            names(fig2a),
            function(x){
                hist(fig2a[[x]],plot=FALSE,breaks=breaks)
            },
            simplify=FALSE
        )
        ylim=range(unlist(sapply(densities,function(x){x$y})),unlist(sapply(hists,function(x){x$density})))
        ylim[2]=ceiling(ylim[2]*10+1)/10

        ## For each GDE method, plot density
        plot.new()
        plot.window(xlim=xlim,ylim=ylim)
        box(bty="l")
        sapply(
            names(fig2a),
            function(x){
                plot(hists[[x]],border=cols[x],col=tp(cols[x]),add=TRUE,freq=FALSE)
            }
        )
        
        ## For each GDE method, plot dashed vertical line at the median
        sapply(
            names(fig2a),
            function(x){
                segments(x0=median(fig2a[[x]]),y0=0,y1=0.95*ylim[2],lty=2,lwd=2,col=cols[x])
                text(x=median(fig2a[[x]]),y=0.98*ylim[2],srt=45,adj=c(0,0.5),labels=signif(median(fig2a[[x]]),2),xpd=NA,col=cols[x])
            }
        )
        
        ## For each GDE method, plot density
        sapply(
            names(fig2a),
            function(x){
                lines(densities[[x]],col=tp(cols[x],"FF"),lwd=0.75)
            }
        )

        ## Add axes, labels, titles
        axis(side=1,at=seq(xlim[1],xlim[2],by=1))
        title(xlab="-log10(P)",line=2,xpd=NA,ylab="Density")
        title(main=tools::toTitleCase(gsub("X_","",i)))
        axis(side=2)

        ## Add color code
        if(i==names(ps)[1]){
            legend(x="topright",lty=1,legend=names(cols),col=cols,bty="n",cex=0.8,inset=c(0,0))
        }
    }
)
dev.off()

####################
## Panels D) and E) data prep
####################

## Load results. The script that generated them is in logistic_regression/figshare/EB/simulations_analysis_clean.R
splatter_dir="./figshare/splatter/"
file=file.path(splatter_dir,"rds","simulations_pre_degs.rds")
## files=list.files(file.path(splatter_dir,"rds/"),full.names=FALSE,include.dirs=FALSE)
## files=files[!grepl("simulations_pre_degs.rds",files)]
files=paste0(c("authors","authors_highES","vanilla"),"_degs.rds")
## Load all simulations' results
simulations=sapply(
    files,
    function(x){
        print(x)
        readRDS(file.path(splatter_dir,"rds",x))
    },
    simplify=FALSE
)

## Extract DEG results and harmonize names. MAST/LR/AOV are statistical methods, TX or GN is about whether the analysis is done at the transcript or gene level, and sidak/multiv means that pvalues are aggregated from transcript to gene levels, univ means univariate analysis
## simulations=lapply(
##     simulations,
##     function(x){
##         degs=x$degs
##         x$degs_formatted=list(
##             mast_tx_sidak=degs$mast_tpm$gn,
##             mast_tx_univ=degs$mast_tpm$tx,
##             mast_gn_univ=degs$mast_unaggregated_gn_tpm$gn,
##             lr_tx_multiv=degs$lr_tpm,
##             lr_tx_univ=degs$lr_unaggregated_tpm,
##             lr_gn_univ=degs$lr_unaggregated_gn_tpm,
##             aov_tx_sidak=degs$aov_tpm$gn,
##             aov_tx_univ=degs$aov_tpm$tx,
##             aov_gn_univ=degs$aov_unaggregated_gn_tpm$gn
##         )
##         x
##     }
## )

simulations = lapply(
    simulations,
    function(x){
        degs = x$degs
        x$degs_formatted = list(
            mast_tx_sidak=degs$mast_tpm$gn,
            mast_tx_univ=degs$mast_tpm$tx,
            mast_gn_univ=degs$mast_unaggregated_gn_tpm$gn,
            lr_tx_multiv=degs$lr_tpm,
            lr_tx_univ=degs$lr_unaggregated_tpm,
            lr_gn_univ=degs$lr_unaggregated_gn_tpm,
            edgeR_tx_sidak=degs$edgeRQLF_CDR$gn,
            edgeR_tx_univ=degs$edgeRQLF_CDR$tx,
            edgeR_gn_univ=degs$edgeRQLF_CDR_gn,
            limma_tx_sidak=degs$limmatrend$gn,
            limma_tx_univ=degs$limmatrend$tx,
            limma_gn_univ=degs$limmatrend_gn,
            DESeq2_tx_sidak=degs$DESeq2$gn,
            DESeq2_tx_univ=degs$DESeq2$tx,
            DESeq2_gn_univ=degs$DESeq2_gn
        )
        x
    }
)

## For all simulations
##     For all pvalue thresholds
##          Compute TPR, FPR, total number of discoveries
all_sim=sapply(
    simulations,
    function(x){
        g=names(x$degs_formatted$lr_tx_multiv) ## Make sure that all vectors are in the same order
        de_genes=sapply(
            split_matrix(as.matrix(rowData(x$splat)),rowData(x$splat)$gene_ensembl),
            function(y){
                any(rowSums(as.matrix(y[,grepl("DEFacGroup",colnames(y))])!=1)>0)
            }
        )
        TP_max=sum(de_genes)
        TN_max=sum(!de_genes)

        ## Reshape data
        df=cbind(
            lr_tx_multiv=x$degs_formatted$lr_tx_multiv[g],
            lr_gn_univ=x$degs_formatted$lr_gn_univ[g],
            mast_tx_sidak=x$degs_formatted$mast_tx_sidak[g],
            mast_gn_univ=x$degs_formatted$mast_gn_univ[g],
            edgeR_tx_sidak=x$degs_formatted$edgeR_tx_sidak[g],
            edgeR_gn_univ=x$degs_formatted$edgeR_gn_univ[g],
            limma_tx_sidak=x$degs_formatted$limma_tx_sidak[g],
            limma_gn_univ=x$degs_formatted$limma_gn_univ[g],
            DESeq2_tx_sidak=x$degs_formatted$DESeq2_tx_sidak[g],
            DESeq2_gn_univ=x$degs_formatted$DESeq2_gn_univ[g],
            DE=de_genes[g]
        )
        df=as.data.table(df)
        ##df=melt(df,measure.vars=c("lr_tx_multiv","mast_tx_sidak","mast_gn_univ"),id.vars="DE")
        df=melt(df,measure.vars=setdiff(colnames(df), "DE"),id.vars="DE")
        colnames(df)=c("DE","method","p")
        df[,simulation:=rep(x$name,nrow(df))]
        df=split(df,as.character(df$method))

        ## For each method, for each p-value threshold, compute TPR and FPR
        df=lapply(df,function(y){
            y[,q:=p.adjust(p,method="fdr")]
            y[,n:=nrow(.SD),by=p]

            y[,TPR:=sum(.SD$DE),by=p]
            y[,FPR:=sum(!.SD$DE),by=p]
            setorder(y,p)
            y=y[!duplicated(p),]
            y=y[,n:=cumsum(y$n)]
            y=y[,TPR:=cumsum(y$TPR)/TP_max]
            y=y[,FPR:=cumsum(y$FPR)/TN_max]
            y
        })
        df=do.call(rbind,df)
        df
    },simplify=FALSE
)
all_sim=do.call(rbind,all_sim)
all_sim[,method:=as.character(all_sim$method)]

## Choose colors and linetypes
methods=unique(all_sim$method)
## colors=c(
##     lr_tx_multiv="#A6CEE3",
##     mast_tx_sidak="#1F78B4",
##     lr_gn_univ="#B2DF8A",
##     mast_gn_univ="#33A02C"
## )
colors = setNames(
    brewer.pal(length(methods), "Paired")[c(7:8,3:4,9:10,1:2,5:6)],
    methods
)

linetype=setNames(rep("solid",length(methods)),methods)

####################
## Panel D): Simulations, sensitivity vs nominal FDR (q-value) threshold
####################

simulations_selected=c("authors","authors_highES")
colours_selected=setNames(brewer.pal(10,"Paired")[7:10],paste(rep(simulations_selected,each=2),rep(c("lr_tx_multiv","mast_tx_sidak"),2)))
colours_selected=c(colours_selected,setNames(brewer.pal(10,"Paired")[c(8,10)],paste(simulations_selected,"mast_gn_univ")))

df_tmp=subset(all_sim,simulation%in%simulations_selected&q<0.05)
df_tmp=do.call(rbind,lapply(split(df_tmp,paste0(df_tmp$simulation,df_tmp$method)),function(x){x[which.max(q),]}))

####################
## New panels C: Show that MAST benefits from aggregation
####################

df3=subset(all_sim, simulation %in% simulations_selected & method %in% c("mast_tx_sidak", "mast_gn_univ"))
df3 = subset(all_sim, simulation %in% "authors_highES")
df3ag=do.call(rbind,lapply(split(df3,paste0(df3$simulation,df3$method)),function(x){
    x=subset(x,q<0.05)
    x[which.max(q),]
}))
groups=with(df3ag,paste(simulation,method))
## linetypes=setNames(ifelse(grepl("mast_gn_univ",groups),"dashed","solid"),groups)
linetypes = setNames(rep("solid",length(methods)),methods)

aucs_1a = df3[, list(auc = auc(TPR, FPR)), by = c("method", "simulation")]
file=file.path("./graphs/table_1a.csv")
write.csv(aucs_1a, file)


p3=ggplot(data=df3,aes(x=FPR,y=TPR,linetype=method,colour=method,group=method)) +
    geom_path() +
    xlab("FPR") +
    ylab("TPR") +
    theme_bw() +
    xlim(0,1) +
    ylim(0,1) +
    geom_abline(intercept=0,slope=1,linetype="dashed") +
    ## ggtitle("Simulations: ROC curve") +
    ##theme(legend.position = "none") +
    scale_colour_manual(values=colors) +
    scale_linetype_manual(values=linetypes) +
    geom_point(data=df3ag,aes(x=FPR,y=TPR,colour=method,group=method))

pdf("./graphs/1A.pdf",height=3.5,width=5)
p3
dev.off()

####################
## Analysis of experimental datasets
####################

## Load results. These are generated in logistic_regression/figshare/EB/shuffling_analysis.R
file=file.path(splatter_dir,"shuffle_rds","neg.rds")
data=readRDS(file)

## We use the analysis where transcripts present in less than 20% of the cells are filtered out
cutoff=as.character(0.2)
## Load data and harmonize
df=lapply(
    names(data),
    function(ds){
        ## df=data.table(
        ##     lr_noCDR_noscrbl=data[[ds]]$degs[[as.character(cutoff)]]$lr_noCDR_noscrbl$gn,
        ##     mast_CDR_noscrbl=data[[ds]]$degs[[as.character(cutoff)]]$mast_CDR_noscrbl$gn,
        ##     lr_noCDR_scrbl=data[[ds]]$degs[[as.character(cutoff)]]$lr_noCDR_scrbl$gn,
        ##     mast_CDR_scrbl=data[[ds]]$degs[[as.character(cutoff)]]$mast_CDR_scrbl$gn
        ## )

        df=data.table(
            lr_noscrbl=data[[ds]]$degs[[cutoff]]$lr_noCDR_noscrbl,
            lr_scrbl=data[[ds]]$degs[[cutoff]]$lr_noCDR_scrbl,
            mast_CDR_noscrbl=data[[ds]]$degs[[cutoff]]$mast_CDR_noscrbl$gn,
            mast_CDR_scrbl=data[[ds]]$degs[[cutoff]]$mast_CDR_scrbl$gn,
            edgeR_noscrbl=data[[ds]]$degs[[cutoff]]$edgeRQLF_CDR_noscrbl$gn,
            edgeR_scrbl=data[[ds]]$degs[[cutoff]]$edgeRQLF_CDR_scrbl$gn,
            limma_noscrbl=data[[ds]]$degs[[cutoff]]$limmatrend_noscrbl$gn,
            limma_scrbl=data[[ds]]$degs[[cutoff]]$limmatrend_scrbl$gn,
            DESeq2_noscrbl=data[[ds]]$degs[[cutoff]]$DESeq2_noscrbl$gn,
            DESeq2_scrbl=data[[ds]]$degs[[cutoff]]$DESeq2_scrbl$gn
        )
        
        ## Compute freq. of discovery
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
            x=rbind(x,data.table(method=x$method[1],p=1,o=max(x$o)+1,q=1,n=max(x$n)))
            x[,dataset:=ds]
            x
        })
        
        df
    }
)

power_1f = lapply(
        df,
        function(x){
            sapply(
                x,
                function(y){
                    suppressWarnings(y[q<0.05, max(n)])
                }
            )
        }
)
names(power_1f) = sapply(
    df,
    function(x){
        x[[1]]$dataset[1]
    }
)
fdr_1g = sapply(
    power_1f,
    function(x){
        x[grepl("_scrbl", names(x))]
    }
)
fdr_1g[!is.finite(fdr_1g)] = 0
power_1f = sapply(
    power_1f,
    function(x){
        x[grepl("_noscrbl", names(x))]
    }
)
file=file.path("./graphs/table_1f.csv")
write.csv(power_1f, file)
file=file.path("./graphs/table_1g.csv")
write.csv(fdr_1g, file)

df=do.call(rbind,lapply(df,function(x){do.call(rbind,x)}))

write.table(power_1f, file 

## Define colors and linetypes
datasets = unique(df$dataset)
methods=apply(expand.grid(c("mast_CDR","lr", "edgeR", "limma", "DESeq2"),c("_scrbl","_noscrbl")),1,paste0,collapse="")
colors = brewer.pal(5, "Set1")
colors = setNames(colors, methods[grepl("_noscrbl", methods)])
## colors=apply(unique(df[,c("dataset","method")]),1,paste,collapse=" ")
## colors=colors[grepl("_noscrbl",colors)]
## colors=setNames(do.call(c,lapply(split(brewer.pal(10,"Paired")[1:6], rep(1:3, each = 2)),function(x){colorRampPalette(x)(5)})), colors)
colors_scrbl=colors
names(colors_scrbl)=sub("noscrbl","scrbl",names(colors_scrbl))

dataset_names = c(
    "Authors' 10X data" = "Ntranos",
    "Petropoulos et al's data" = "Petropoulos",
    "GSE64016" = "Leng",
)

## Plot
p3 = lapply(
    datasets,
    function(dataset){
        envir = environment()
        ggplot(data=df[grepl("_noscrbl",method) & dataset == get("dataset", envir = envir),],aes(x=q, y=n, colour=method, group=paste(dataset,method))) +
            geom_vline(xintercept=0.05,linetype="dashed") +
            geom_path(size=0.5) +
            theme_bw() +
            ggtitle(dataset) +
            xlab("Nominal FDR") +
            ylab("Fraction of discoveries") +
            scale_colour_manual(values=colors) +
            theme(legend.position="none")
    }
)

p4 = lapply(
    datasets,
    function(dataset){
        envir = environment()
        ggplot(data=df[grepl("_scrbl",method) & dataset == get("dataset", envir = envir),],aes(x=q, y=n, colour=method, group=paste(dataset,method))) +
            geom_vline(xintercept=0.05,linetype="dashed") +
            geom_path(size=0.5) +
            theme_bw() +
            ggtitle(dataset) +
            xlab("Nominal FDR") +
            ylab("Fraction of discoveries") +
            scale_colour_manual(values=colors_scrbl) +
            theme(legend.position="none")
    }
)

file="./graphs/Panel_ROCs.pdf"
pdf(file=file,width=17/5*3,height=3.5*2)
ggarrange(plotlist = c(p3, p4), nrow = 2, ncol = 3)
dev.off()

pdf(file=sub(".pdf","_legend.pdf", file, fixed = TRUE), ,width=2,height=2)
par(mar = c(0, 0, 0, 0))
plot.new()
ys = seq(0.8, 0.2, length.out = 5)
text(x = 0.2, y = ys, labels = sub("_noscrbl", "", names(colors)), pos = 4)
segments(x0 = 0.18, x1 = 0.1, y0 = ys, col = colors, lwd = 4)
dev.off()

####################
## Correlation of transcripts in simulations vs real data
####################

data_correlation=simulations[paste0(c("authors","authors_highES","vanilla"),"_degs.rds")]
data_correlation=lapply(
    data_correlation,
    function(x){
        xp=log2(1+assays(x$splat)$count)
        expressed=rowSums(xp)>0
        xp=xp[expressed,]
        t2g=rowData(x$splat)$gene_ensembl[expressed]
        x=list(xp=xp,t2g=t2g,name=x$name)
    }
)
data_correlation=c(
    data_correlation,
    lapply(
        data,
        function(x){
            expressed=Matrix:::rowSums(x$xp)>0
            x$xp=x$xp[expressed,]
            x$t2g=as.character(x$t2g)[expressed]
            x[c("xp","t2g","name")]
        }
    )
)
data_correlation=lapply(
    data_correlation,
    function(x){
        print(x$name)
        xp=x$xp
        t2g=as.character(x$t2g)
        xp_split=split_matrix(as.matrix(xp),t2g)
        xp_split=xp_split[sapply(xp_split,length)>1]
        xp_split=lapply(xp_split,t)
        cors_intra=lapply(
            xp_split,
            cor
        )
        cors_intra=lapply(cors_intra,function(x){x[upper.tri(x)]})

        xp=as.matrix(t(xp))
        spl=sample(1:ncol(xp),min(ncol(xp),2000))
        cors_inter=cor(xp[,spl])
        cors_inter=cors_inter[upper.tri(cors_inter)]
        x=c(x,list(cors_inter=cors_inter,cors_intra=cors_intra))
        x
    }
)

n=ceiling(length(data_correlation)/2)
png(file.path("./graphs/","Sup Fig 4_sqrt.png"),res=300,height=n*1000,width=2000)
layout(matrix(nrow=n,ncol=2,data=1:(3*n),byrow=FALSE))
lapply(
    data_correlation,
    function(x){
        densities_intra=density(unlist(x$cors_intra),from=-1,to=1,n=512)
        densities_inter=density(x$cors_inter,from=-1,to=1,n=512)
        densities_intra=densities_intra[c("x","y")]
        densities_inter=densities_inter[c("x","y")]
        densities_intra$x=c(-1,densities_intra$x,1)
        densities_inter$x=c(-1,densities_inter$x,1)
        densities_intra$y=c(0,densities_intra$y,0)
        densities_inter$y=c(0,densities_inter$y,0)
        ## densities_intra$y=log10(1+densities_intra$y)
        ## densities_inter$y=log10(1+densities_inter$y)
        densities_intra$y=sqrt(densities_intra$y)
        densities_inter$y=sqrt(densities_inter$y)
        ylim=c(0,max(c(densities_intra$y,densities_inter$y)))
        plot.new()
        plot.window(ylim=ylim,xlim=c(-1,1))
        axis(side=1)
        axis(side=2)
        box(bty="l")
        polygon(densities_intra,border="#FF0000CC",col="#FF000044")
        polygon(densities_inter,border="#0000FFCC",col="#0000FF44")
        title(xlab="Pearson's r",ylab="",ylim=ylim)
        title(main=x$name,cex.main=par()$cex.main*0.8,font.main=3,line=0)
        if(x$name==data_correlation[[1]]$name){
            legend(
                x="topright",
                lty=1,
                col=c("#00000000","red","blue"),
                legend=c("Set of transcripts:","of the same gene","randomly selected"),
                bty="n",
                cex=0.9
            )
            title(main="Simulated datasets",line=3)
        }
        if(x$name%in%sapply(data_correlation,function(y){y$name})[4]){
            title(main="Real datasets",line=3)            
        }
        title(ylab="sqrt(Density)")
    }
)
dev.off()

####################
## Simulations modeled from new datasets
####################

##
## Reformat DEGs (for consistency)
##

simulations=readRDS("../figshare/splatEstimate/rds/simulations.rds")
names_map=c(
    "TM_thymus"="Tabula Muris - (10x v2)",
    "10xv2PBMC4k"="PBMC4k (10X v2)",
    "10xv3PBMC1k"="PBMC1k (10X v3)",
    "EMTAB3929"="Petropoulos (Smart-Seq2)"
)

plots=lapply(
    simulations[names(names_map)],
    function(x){
        
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
        methods = c(
            "mast_tx_sidak",
            "lr_tx_multiv",
            "edgercdr_tx_sidak",
            "limma_tx_sidak",
            "deseq2_tx_sidak"
        )
        methods = rev(methods)
        
        fdrs=sapply(
            methods
           ,function(y){
               if(grepl("ENSG",names(x$degs_formatted[[y]])[1])){
                   true_DE=de_genes
               } else if (grepl("ENST",names(x$degs_formatted[[y]])[1])){
                   true_DE=de_transcripts
               }
               calculate_fdr(true_DE=true_DE[names(x$degs_formatted[[y]])],pvalues=x$degs_formatted[[y]],title=y)
           },simplify=FALSE)

        fdrs <- lapply(
            1:length(fdrs),
            function(y){
                data.frame(
                    sensitivity=fdrs[[y]]$sen,
                    actual_FDR=fdrs[[y]]$fdr,
                    FalsePositiveRate = fdrs[[y]]$FalsePositiveRate,
                    nominal_FDR=p.adjust(fdrs[[y]]$pvalues,method="fdr"),
                    Method = names(fdrs)[y])
            }
        )
        fdrs <- do.call(rbind, fdrs)

        aucs = sapply(
            split(fdrs, fdrs$Method),
            function(df){
                ## i = 2:nrow(df) - 1
                ## heights = (df[i, "sensitivity"]+df[i+1, "sensitivity"])/2
                ## widths = df[i+1, "FalsePositiveRate"]-df[i, "FalsePositiveRate"]
                ## sum(heights*widths)
                auc(TPR = df[, "sensitivity"], df[, "FalsePositiveRate"])
            }
        )
        aucs_FPR10 = sapply(
            split(fdrs, fdrs$Method),
            function(df){
                df = df[df$FalsePositiveRate<0.1, ]
                auc(TPR = df[, "sensitivity"], df[, "FalsePositiveRate"])*10
            }
        )
        
        ## colors=setNames(makeTransparent(rep(brewer.pal(length(methods)/2,"Set1"),each=2),123),unique(fdrs$Method))
        ## linetype=setNames(rep(c("solid","dashed","dotted")[1:2],length(methods)/2),unique(fdrs$Method))
        colors=setNames(makeTransparent(rev(rep(brewer.pal(length(methods),"Set1"),each=1)),255),unique(fdrs$Method))
        linetype=setNames(rep(c("solid","dashed","dotted")[1],length(methods)),unique(fdrs$Method))

        p2 <- ggplot(
            data = fdrs,
            aes(x = FalsePositiveRate , y = sensitivity, colour=Method, linetype=Method)
        ) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            geom_path(size = .5) +
            ## geom_abline(intercept=0,slope=1,linetype="solid",col="gray") +
            theme_bw(12) +
            coord_fixed() +
            scale_x_sqrt() +
            ggtitle(names_map[x$name]) +
            theme(plot.title = element_text(hjust = 0.5)) +
            xlab("FPR") +
            ylab("TPR")

        ## zoomtheme <- theme(legend.position="none", axis.line=element_blank(),axis.text.x=element_blank(),
        ##     axis.text.y=element_blank(),axis.ticks=element_blank(),
        ##     axis.title.x=element_blank(),axis.title.y=element_blank(),
        ##     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        ##     panel.background = element_rect(color='red', fill="white"),
        ##     plot.margin = unit(c(0,0,-6,-6),"mm"))

        ## p2_zm = p2 + zoomtheme
        ## vp = viewport(1, 0, adjust = c())
        list(plot = p2, aucs = aucs, aucs_FPR10 = aucs_FPR10)
    }
)

file=file.path("./graphs/1E_new_simulations.pdf")
library(RSvgDevice)
pdf(file,height=3.5,width=8*4/3)
ggarrange(plotlist=lapply(plots, function(x){x[["plot"]]}),ncol=4,common.legend=TRUE)
dev.off()

file=file.path("./graphs/table_1e.csv")
aucs_1e = do.call(
    rbind,
    lapply(
        plots,
        function(x){
            df1 = as.data.frame(as.list(x[["aucs"]]))
            df2 = as.data.frame(as.list(x[["aucs_FPR10"]]))
            colnames(df1) = paste("AUC", colnames(df1))
            colnames(df2) = paste("AUC_FPR10", colnames(df2))
            cbind(
                df1,
                df2
            )
        }
    )
)
rownames(aucs_1e) = names_map[rownames(aucs_1e)]
write.csv(signif(aucs_1e, 2), file = file)
