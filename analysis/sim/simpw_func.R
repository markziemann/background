#library("tidyverse")
library("parallel")
library("edgeR")
library("DESeq2")
library("limma")
library("stringi")
library("mitch")
library("fgsea")
library("clusterProfiler")

########################################
# get some counts
########################################
countData<-function() {
# a is orig expression data
a<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/ERR2539161.se.tsv")

# GI is the gene information. It has gene name information for mapping accessions to gene symbols in GMT files
gi<-read.table("https://raw.githubusercontent.com/markziemann/simde/master/start_data/ERR2539161/GeneInfo.tsv",header=T,row.names=1)

# merge gene names
aa<-merge(a,gi,by=0)
aa<-aa[,-c(4:7)]
aa<-aggregate(. ~ GeneSymbol,aa,function(x) sum(as.numeric(as.character(x))))
aa$Row.names=NULL
rownames(aa)<-aa$GeneSymbol
aa$GeneSymbol=NULL
a<-aa[which(aa$ERR2539161>=10),,drop=F]
a
}

########################################
# generate some gene sets
########################################
randomGeneSets<-function(a,setsize,nsets){
gsets<-sapply( rep(setsize,nsets) , function(x) {list(as.character(sample(rownames(a),x))) } )
names(gsets)<-stri_rand_strings(length(gsets), 15, pattern = "[A-Za-z]")
gsets
}

########################################
# simulate some gene expression data
########################################
simrna<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets) {

# N_REPS=5 ; SUM_COUNT=10000000 ; VARIANCE=0.2 ; FRAC_DE=0.05 ; FC=1 ; GMT="ReactomePathways.gmt"

library("edgeR")

df = NULL
for (k in paste0("data",1:(N_REPS*2)))  {
        b<-thinCounts(a,target.size=SUM_COUNT)
        colnames(b)=k
        df = cbind(df,b)
     }

# now need to only include gsets with 10 members in the 
gsets_sub<-which(unlist( lapply(gsets,function(x) { length(which(rownames(a) %in% as.character(unlist(x)))) >10 }  ) ) )
gsets<-gsets[which(names(gsets) %in% names(gsets_sub))]

#Number of differential gene sets
NDIF=round(length(gsets)*FRAC_DE)

if (VARIANCE>0) {
  #create some random values centred around 1 with some% error
  rand<-matrix(log2(rnorm(nrow(a)*N_REPS*2 , 2, VARIANCE)),ncol=N_REPS*2)
  #incorporate the noise
  df<-round(df*rand)
  #set any negative counts to zero
  df<-apply(df, 2, function(x) {ifelse(x < 0, 0, x)})
} 

if (NDIF>0) {
  message("prep fold changes")
  #Make even
  if ( NDIF%%2==1 ) { print("odd") ; NDIF=NDIF-1 }

  # sample some pathways to fiddle with
  DE_LIST<-sample(gsets , NDIF)

  # divide the list in 2 with half up and half down
  UP_LIST=sample(DE_LIST , NDIF/2)
  DN_LIST<-DE_LIST[!(DE_LIST %in% UP_LIST)]

  # now find a list of genes inside the pathways
  UP_DE<-unique(unlist(unname(UP_LIST)))

  # select the ones that are also in the profile
  UP_DE<-UP_DE[which(UP_DE %in% row.names(df))]

  # same for down genes
  DN_DE<-unique(unlist(unname(DN_LIST)))
  DN_DE<-DN_DE[which(DN_DE %in% row.names(df))]

  ITX<-intersect(UP_DE,DN_DE)
  # need to eliminate the overlapping ones for simplicity
  UP_DE<-setdiff(UP_DE,ITX)
  DN_DE<-setdiff(DN_DE,ITX)

  #reformat as df and add fold change
  UP_DE<-as.data.frame(UP_DE)
  UP_DE$V1<-2^FC
  colnames(UP_DE)=c("Gene","FC")
  DN_DE<-as.data.frame(DN_DE)
  DN_DE$V1<-2^-FC
  colnames(DN_DE)=c("Gene","FC")
  ALL_DE<-rbind(DN_DE,UP_DE)
  #Go back to list for downstream work
  UP_DE<-UP_DE$Gene
  DN_DE<-DN_DE$Gene
  NON_DE<-as.data.frame(setdiff(rownames(df),ALL_DE$Gene))
  colnames(NON_DE)="Gene"
  NON_DE$FC=1
  ALL_DE<-rbind(ALL_DE,NON_DE)
  ALL_DE<-ALL_DE[ order(as.vector(ALL_DE$Gene)) , ]
  message("incorporate changes")

  df <- df[ order(row.names(df)), ]
  df2<-cbind(df,ALL_DE)
  df2$Gene=NULL
} else {
  df2<-as.data.frame( df )
  df2$FC<- 1
  UP_DE=NULL
  DN_DE=NULL
  UP_LIST=NULL
  DN_LIST=NULL
}
ODD_COLS=(1:(ncol(df2)-1))[c(TRUE,FALSE)]
EVEN_COLS=(1:(ncol(df2)-1))[c(FALSE,TRUE)]
controls<-df2[,ODD_COLS]
colnames(controls)=paste0( "ctrl_" ,1:ncol(controls) )
treatments<-round(df2[,EVEN_COLS]*df2$FC)
colnames(treatments)=paste0( "trt_" ,1:ncol(treatments) )
x<-cbind(controls,treatments)
rownames(x)=rownames(df2)
#filter out genes that are not expressed
x<- x[which(rowSums(x)/ncol(x)>10),]
UP_DE<-intersect(UP_DE,rownames(x))
DN_DE<-intersect(DN_DE,rownames(x))
xx <- list("x" = x,"UP_DE"=UP_DE,"DN_DE"=DN_DE,"UP_LIST"=UP_LIST,"DN_LIST"=DN_LIST)
xx
}

#################################################
# a parallel repeat function
##################################################
#Thanks Gray Calhoun gcalhoun@iastate.edu for the following function
RepParallel <- function(n, expr, simplify = "array",...) {
      answer <-
        mclapply(integer(n), eval.parent(substitute(function(...) expr)),...)
      if (!identical(simplify, FALSE) && length(answer)) 
        return(simplify2array(answer, higher = (simplify == "array")))
      else return(answer)
    }
# RepParallel usage
#xxx<-RepParallel(10,simrna(a,5,10000000,0.2,20), simplify=F, mc.cores = detectCores() )


#################################################
# define DESeq2 function
##################################################
deseq2<-function(x) {
  library("DESeq2")
  label = "simulate"
  y <- x[[1]]
  samplesheet <- as.data.frame(colnames(y))
  colnames(samplesheet) = "sample"
  samplesheet$trt <- factor(as.numeric(grepl("trt",colnames(y))))
  dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
  res <- DESeq(dds)
  z <- DESeq2::results(res)
  x[[6]] <- as.data.frame(z[order(z$pvalue),])
  up <- rownames(subset(z,log2FoldChange>0 & padj<0.05))
  if (length(up)<250) { up <- rownames(head(subset(z,log2FoldChange>0 ),250)) }
  dn <- rownames(subset(z,log2FoldChange<0 & padj<0.05))
  if (length(dn)<250) { dn <- rownames(head(subset(z,log2FoldChange<0 ),250)) }
  x[[7]] <- up
  x[[8]] <- dn
  x
}


#################################################
# define mitch function
##################################################
run_mitch <- function(y,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=16) {
  library("mitch")
  dge <- sapply(y,"[",6)
  names(dge) <- paste0("x",1:length(dge),sep="")

  mres <- mclapply(1:length(dge),function(d) {
    w <- mitch_import(dge[[d]], DGE_FUNC)
    mres <- mitch_calc(w,gsets,priority="significance",cores=1,minsetsize=5)
  } , mc.cores=cores)

  obs_up <- lapply(1:length(mres), function(d) {
        up <- subset(mres[[d]]$enrichment_result, p.adjustANOVA<0.05 & s.dist>0)[,1]
        y[[d]][[19]] <- up
        up
  } )

  obs_dn <- lapply(1:length(mres), function(d) {
        dn <- subset(mres[[d]]$enrichment_result, p.adjustANOVA<0.05 & s.dist<0)[,1]
        y[[d]][[20]] <- dn
        dn
  } )

  gt_up <- sapply(y,"[",4)
  gt_up <- lapply( gt_up , names)
  gt_dn <- sapply(y,"[",5)
  gt_dn <- lapply( gt_dn , names)

  true_pos_up <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
  true_pos_dn <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
  false_pos_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
  false_pos_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
  false_neg_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
  false_neg_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

  true_pos <- mean(true_pos_up+true_pos_dn)
  false_pos <- mean(false_pos_up+false_pos_dn)
  false_neg <- mean(false_neg_up+false_neg_dn)
  nrows <- as.numeric(lapply( sapply(y,"[",1 ), nrow))
  true_neg <- mean(nrows-(true_pos+false_pos+false_neg))

  p <- true_pos/(true_pos+false_pos)
  r <- true_pos/(true_pos+false_neg)
  f <- 2*p*r/(p+r)

  attr(y,'mitch') <- data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
  y
}


##################################
# clusterprofiler default function
##################################
# Note that clusterprofiler requires different gene set format
run_clusterprofiler_default <- function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=16){

gsets2 <- stack(gsets)[,c(2,1)]
colnames(gsets2) <- c("term","gene")

dge <- sapply(x,"[",6)
up <- sapply(x,"[",7)
dn <- sapply(x,"[",8)
bgs <- lapply(dge, function(x) { rownames(x) } )

l_ups <- sapply(ups,length)
l_dns <- sapply(dns,length)
geneset_sizes <- sapply( gsets2 , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

obs_up <- mclapply(1:length(dge), function(d) {
  if ( length(which(ups[[d]] %in% gsets2$gene)) > 0 ) {
    
    ora_up <- as.data.frame(enricher(gene = ups[[d]] ,
      universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
      pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

    ora_up$geneID <- NULL
    ora_up <- subset(ora_up,p.adjust<0.05 )
    ora_ups <- rownames(ora_up)
    ora_ups
  } else {
    ora_ups = NULL
  }
  x[[d]][[9]] <- ora_ups
  ora_ups
} , mc.cores=cores )

obs_dn <- mclapply(1:length(dge), function(d) {
  if ( length(which(dns[[d]] %in% gsets2$gene)) > 0 ) {

    ora_dn <- as.data.frame(enricher(gene = dns[[d]] ,
      universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
      pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

    ora_dn$geneID <- NULL
    ora_dn <- subset(ora_dn,p.adjust<0.05 )
    ora_dns <- rownames(ora_dn)
  } else {
    ora_dns = NULL
  }
  x[[d]][[10]] <- ora_dns
  ora_dns
} , mc.cores=cores )

#ground truth comparison
gt_up <- sapply(x,"[",4)
gt_up <- lapply( gt_up , names)
gt_dn <- sapply(x,"[",5)
gt_dn <- lapply( gt_dn , names)

true_pos_up <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos <- mean(true_pos_up+true_pos_dn)
false_pos <- mean(false_pos_up+false_pos_dn)
false_neg <- mean(false_neg_up+false_neg_dn)
nrows <- as.numeric(lapply( sapply(x,"[",1 ), nrow))
true_neg <- mean(nrows-(true_pos+false_pos+false_neg))

p <- true_pos/(true_pos+false_pos)
r <- true_pos/(true_pos+false_neg)
f <- 2*p*r/(p+r)

attr(x,'clusterprofiler') <- data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# FORA function
##################################
run_fora <- function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=16){

dge <- sapply(x,"[",6)

bg <- lapply( dge, function(d) {
  rownames(d)
} )

up <- sapply(x,"[",7)

dn <- sapply(x,"[",8)

obs_up <- mclapply( 1:length(dge) , function(d) {
  subset(as.data.frame(fora(pathways=gsets, genes=up[[d]],universe=bg[[d]],minSize = 2)),padj<0.05)$pathway
}, mc.cores=cores)

obs_dn <- mclapply( 1:length(dge) , function(d) {
  subset(as.data.frame(fora(pathways=gsets, genes=dn[[d]],universe=bg[[d]],minSize = 2)),padj<0.05)$pathway
}, mc.cores=cores)

for (d in 1:length(dge)) {
  x[[d]][[15]] <- obs_up[[d]]
  x[[d]][[16]] <- obs_dn[[d]]
}

gt_up <- sapply(x,"[",4)
gt_up <- lapply( gt_up , names)
gt_dn <- sapply(x,"[",5)
gt_dn <- lapply( gt_dn , names)

true_pos_up <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos <- mean(true_pos_up+true_pos_dn)
false_pos <- mean(false_pos_up+false_pos_dn)
false_neg <- mean(false_neg_up+false_neg_dn)
nrows <- as.numeric(lapply( sapply(x,"[",1 ), nrow))
true_neg <- mean(nrows-(true_pos+false_pos+false_neg))

p <- true_pos/(true_pos+false_pos)
r <- true_pos/(true_pos+false_neg)
f <- 2*p*r/(p+r)

attr(x,'fora') <- data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# FGSEA function
##################################
run_fgsea <- function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=16 ){

dge <- sapply(x,"[",6)

xx <- mclapply( dge , function(x) { 
 s <- x$stat
 names(s) <- rownames(x)
 p <- as.data.frame(fgsea(pathways=gsets, stats=s ,minSize = 5 ))
 p
} , mc.cores=cores )

obs_up <- lapply(xx, function(x) { subset(x,padj<0.05 & ES>0)[,1] } )
obs_dn <- lapply(xx, function(x) { subset(x,padj<0.05 & ES<0)[,1] } )

for (d in 1:length(dge)) {
  x[[d]][[17]] <- obs_up[[d]]
  x[[d]][[18]] <- obs_dn[[d]]
}

gt_up <- sapply(x,"[",4)
gt_up <- lapply( gt_up , names)
gt_dn <- sapply(x,"[",5)
gt_dn <- lapply( gt_dn , names)

true_pos_up <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn <- as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn <- as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos <- mean(true_pos_up+true_pos_dn)
false_pos <- mean(false_pos_up+false_pos_dn)
false_neg <- mean(false_neg_up+false_neg_dn)
nrows <- as.numeric(lapply( sapply(x,"[",1 ), nrow))
true_neg <- mean(nrows-(true_pos+false_pos+false_neg))

p <- true_pos/(true_pos+false_pos)
r <- true_pos/(true_pos+false_neg)
f <- 2*p*r/(p+r)

attr(x,'fgsea') <- data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# geneSetTest function
##################################
run_gst<-function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge<-sapply(x,"[",6)

gst_func<-function(gene_names,gset,stats) {
 i<-which(gene_names %in% gset )
 y=geneSetTest(i,stats,alternative = "either",ranks.only=TRUE)
 y=p.adjust(y,method="BH")
 es=mean(rank(s)[i])-mean(rank(stats))
 c(y,es)
}

mygst=NULL
mygst=list()
for (d in 1:length(dge)) {
  s<-dge[[d]]$stat
  mygst[[d]]<-mclapply(gsets, gst_func , gene_names=rownames(dge[[d]]) , stats=s, mc.cores=8)
  mygst[[d]]<-t(as.data.frame(mygst[[d]]))
  colnames(mygst[[d]])<-c("padj","es")
  x[[d]][[15]]<-names(which(mygst[[d]][,1]<0.05 & mygst[[d]][,2]>0))
  x[[d]][[16]]<-names(which(mygst[[d]][,1]<0.05 & mygst[[d]][,2]<0))
}

obs_up<-sapply(x,"[",9)
obs_dn<-sapply(x,"[",10)

gt_up<-sapply(x,"[",4)
gt_up<-lapply( gt_up , names)
gt_dn<-sapply(x,"[",5)
gt_dn<-lapply( gt_dn , names)

true_pos_up<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_up ,  gt_up ))
true_pos_dn<-as.numeric(mapply( function(x,y) length(intersect(x,y)) , obs_dn ,  gt_dn ))
false_pos_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_up ,  gt_up ))
false_pos_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , obs_dn , gt_dn ))
false_neg_up<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_up ,  obs_up ))
false_neg_dn<-as.numeric(mapply( function(x,y) length(setdiff(x,y)) , gt_dn ,  obs_dn ))

true_pos<-mean(true_pos_up+true_pos_dn)
false_pos<-mean(false_pos_up+false_pos_dn)
false_neg<-mean(false_neg_up+false_neg_dn)
nrows<-as.numeric(lapply( sapply(x,"[",1 ), nrow))
true_neg<-mean(nrows-(true_pos+false_pos+false_neg))

p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

attr(x,'gst') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}


##################################
# aggregate function
##################################
agg_dge <- function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets, cores = 32) {

#TEST# N_REPS=5 ; SUM_COUNT=30000000 ; VARIANCE=0.45 ; FRAC_DE=0.05 ; FC=1 ; SIMS=8 ; DGE_FUNC="deseq2" ; gsets=gsets

xxx <- RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets), simplify=F, mc.cores = cores )

# run deseq2
xxx <- mclapply(xxx , DGE_FUNC , mc.cores = cores )

# run clusterprofiler default pos 9,10
xxx <- run_clusterprofiler_default(x=xxx,DGE_FUNC,gsets,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=32)

# run clusterprofiler default pos 11,12
#xxx <- run_clusterprofiler_fixed(x=xxx,DGE_FUNC,gsets,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run clusterprofiler nobg pos 13,14
#xxx <- run_clusterprofiler_nobg(x=xxx,DGE_FUNC,gsets,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run fora pos 15,16
xxx <- run_fora(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=32)

# run fgsea pos 17,18
xxx <- run_fgsea(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=32)

# run mitch pos 19,20
#xxx <- run_mitch(xxx,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,cores=32)


# return the result
g=list()
for (f in 1:length(attributes(xxx))) {
  PWAY_FUNC<-names(attributes(xxx)[f])
  PWAY_FUNC<-as.data.frame(PWAY_FUNC)
  g[[f]]<-cbind(unname(attributes(xxx)[f]),PWAY_FUNC)
}

g
}
# x<-agg_dge(a,10,40000000,0.4,0.2,1,10,"deseq2",gsets) 

