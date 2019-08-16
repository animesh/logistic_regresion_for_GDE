library(patchwork)
library(gplots)
library(data.table)
library(ebmisc)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(gplots)
library(tools)
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
load("../NYMP_2018/10x_example-logR/rdata/Fig2.RData")
load("../NYMP_2018/10x_example-logR/rdata/Fig2_allcells.RData")

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
splatter_dir="../figshare/splatter/"
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
        df=cbind(lr_tx_multiv=x$degs_formatted$lr_tx_multiv[g],mast_tx_sidak=x$degs_formatted$mast_tx_sidak[g],mast_gn_univ=x$degs_formatted$mast_gn_univ[g],DE=de_genes[g])
        df=as.data.table(df)
        df=melt(df,measure.vars=c("lr_tx_multiv","mast_tx_sidak","mast_gn_univ"),id.vars="DE")
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

## We remove the two simulations where we used the "dropout.type="group"" argument in splatter (see https://rdrr.io/github/Oshlack/splatter/man/SplatParams.html). The reason is that this simulates differential dropout across celltypes (so for each gene, cells of celltype A will have e.g. less dropout than cells of celltype B). This is not considered differential expression, although it arguably is. It therefore impacts TPR and FPR in a way we disagree with.
all_sim=subset(all_sim,!grepl("dropout",all_sim$simulation))

## Choose colors and linetypes
methods=unique(all_sim$method)
colors=c(
    lr_tx_multiv="#A6CEE3",
    mast_tx_sidak="#1F78B4",
    lr_gn_univ="#B2DF8A",
    mast_gn_univ="#33A02C"
)
linetype=setNames(rep("solid",length(methods)),methods)

## For each method, interpolate TPR(p) and FPR(p) to a common set of p-values, so that we can plot method A versus method B
interp_points=seq(0,1,by=0.0001)
all_sim_interp=split(all_sim,paste0(all_sim$simulation,all_sim$method))
all_sim_interp=lapply(all_sim_interp,function(x){
    res=data.table(p=interp_points)
    res[,TPR:=approx(x=x$p,y=x$TPR,xout=interp_points)$y] ## Where the interpolation happens
    res[,FPR:=approx(x=x$p,y=x$FPR,xout=interp_points)$y]
    res[p==0,c("TPR","FPR"):=list(0,0)] ## Add points (0,0) and (1,1), which are mathematically guaranteed. Avoid edge cases in the approx() function (interpolating outside of data range)
    res[p==1,c("TPR","FPR"):=list(1,1)]
    res[,method:=rep(x$method[1],nrow(res))]
    res[,simulation:=rep(x$simulation[1],nrow(res))]
    res
})

## Set-up MAST-Sidak vs mLR direct comparison on interpolated data
m1="mast_tx_sidak"
m2="lr_tx_multiv"
all_sim_interp=do.call(rbind,all_sim_interp)
all_sim_interp=split(all_sim_interp,all_sim_interp$simulation)
all_sim_interp=lapply(all_sim_interp,function(x){
    d1=x[method==m1,]
    d2=x[method==m2,]
    colnames(d1)[!colnames(d1)%in%c("p","simulation")]=paste0(colnames(d1)[!colnames(d1)%in%c("p","simulation")],".",m1)
    colnames(d2)[!colnames(d2)%in%c("p","simulation")]=paste0(colnames(d2)[!colnames(d2)%in%c("p","simulation")],".",m2)
    d1=d1[,grep("method.",colnames(d1),value=TRUE):=NULL]
    d2=d2[,grep("method.",colnames(d2),value=TRUE):=NULL]
    merge(d1,d2,by=c("p","simulation"))
})

all_sim_interp=do.call(rbind,all_sim_interp)

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
df3ag=do.call(rbind,lapply(split(df3,paste0(df3$simulation,df3$method)),function(x){
    x=subset(x,q<0.05)
    x[which.max(q),]
}))
groups=with(df3ag,paste(simulation,method))
linetypes=setNames(ifelse(grepl("mast_gn_univ",groups),"dashed","solid"),groups)

p3=ggplot(data=df3,aes(x=FPR,y=TPR,linetype=paste(simulation,method),colour=paste(simulation,method),group=paste(simulation,method))) +
    geom_path() +
    xlab("FPR") +
    ylab("TPR") +
    theme_bw() +
    xlim(0,1) +
    ylim(0,1) +
    geom_abline(intercept=0,slope=1,linetype="dashed") +
    ## ggtitle("Simulations: ROC curve") +
    theme(legend.position = "none") +
    scale_colour_manual(values=colours_selected) +
    scale_linetype_manual(values=linetypes) +
    geom_point(data=df3ag,aes(x=FPR,y=TPR,colour=paste(simulation,method),group=paste(simulation,method)))

pdf("./graphs/1A.pdf",height=3.5,width=3.5)
p3
dev.off()

####################
## Analysis of experimental datasets
####################

## Load results. These are generated in logistic_regression/figshare/EB/shuffling_analysis.R
file=file.path(splatter_dir,"shuffle_rds","neg.rds")
data=readRDS(file)

## We use the analysis where transcripts present in less than 20% of the cells are filtered out
cutoff=0.2
## Load data and harmonize
df=lapply(
    names(data),
    function(ds){
        df=data.table(
            lr_noCDR_noscrbl=data[[ds]]$degs[[as.character(cutoff)]]$lr_noCDR_noscrbl$gn,
            mast_CDR_noscrbl=data[[ds]]$degs[[as.character(cutoff)]]$mast_CDR_noscrbl$gn,
            lr_noCDR_scrbl=data[[ds]]$degs[[as.character(cutoff)]]$lr_noCDR_scrbl$gn,
            mast_CDR_scrbl=data[[ds]]$degs[[as.character(cutoff)]]$mast_CDR_scrbl$gn
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
df=do.call(rbind,lapply(df,function(x){do.call(rbind,x)}))

## Define colors and linetypes
methods=apply(expand.grid(c("mast_CDR","lr_noCDR"),c("_scrbl","_noscrbl")),1,paste0,collapse="")
colors=apply(unique(df[,c("dataset","method")]),1,paste,collapse=" ")
colors=colors[grepl("_noscrbl",colors)]
colors=setNames(brewer.pal(10,"Paired")[1:6],colors)
colors_scrbl=colors
names(colors_scrbl)=sub("noscrbl","scrbl",names(colors_scrbl))

## Plot
p3=ggplot(data=df[grepl("_noscrbl",method),],aes(x=q, y=n, colour=paste(dataset,method),group=paste(dataset,method))) +
    geom_vline(xintercept=0.05,linetype="dashed") +
    geom_path(size=0.5) +
    theme_bw() +
    ggtitle("Real datasets: power analysis") +
    xlab("nominal FDR") +
    ylab("Fraction of discoveries") +
    scale_colour_manual(values=colors) +
    theme(legend.position="none")

p4=ggplot(data=df[grepl("_scrbl",method),],aes(x=q, y=n, colour=paste(dataset,method),group=paste(dataset,method))) +
    geom_vline(xintercept=0.05,linetype="dashed") +
    geom_path(size=0.5) +
    theme_bw() +
    ggtitle("Real datasets, scrambled response: power analysis") +
    xlab("nominal FDR") +
    ylab("Fraction of discoveries") +
    scale_colour_manual(values=colors_scrbl)

file="./graphs/Panel_ROCs.pdf"
pdf(file=file,width=17/5*3,height=3.5)
p3+p4 + plot_layout(nrow=1)
dev.off()

####################
## Correlation of transcripts in simulations vs real data
####################

## Parse simulated and experimental datasets
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

## Compute within-gene and random transcript correlations
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

## Plot
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
    "10xv2PBMC4k"="PBMC4k (10X v2)",
    "10xv3PBMC1k"="PBMC1k (10X v3)",
    "TM_thymus"="Tabula Muris - Thymus"
)

plots=lapply(
    simulations,
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

        ## Compute TPR, FPR...
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

        
        ## Set up ROC curves
        colors=setNames(makeTransparent(rep(brewer.pal(3,"Set1"),each=2),123),unique(fdrs$Method))
        linetype=setNames(rep(c("solid","dashed","dotted")[1:2],3),unique(fdrs$Method))
        p2 <- ggplot(data = fdrs, aes(x = FalsePositiveRate , y = sensitivity, colour=Method, linetype=Method)) +
            scale_colour_manual(values=colors) +
            scale_linetype_manual(values=linetype) +
            ## scale_x_sqrt() +
            geom_path() +
            geom_abline(intercept=0,slope=1,linetype="solid",col="gray") +
            theme_bw(12) +
            coord_fixed() +
            ggtitle(names_map[x$name]) +
            theme(plot.title = element_text(hjust = 0.5)) +
            xlab("FPR") +
            ylab("TPR")
        p2        
    }
)

## Plot ROC curves
file=file.path("./graphs/1E_new_simulations.svg")
library(RSvgDevice)
devSVG(file,height=3.5,width=8)
ggarrange(plotlist=plots,ncol=3,common.legend=TRUE)
dev.off()
