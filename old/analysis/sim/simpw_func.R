library("tidyverse")
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
label="simulate"
y<-x[[1]]
samplesheet<-as.data.frame(colnames(y))
colnames(samplesheet)="sample"
samplesheet$trt<-factor(as.numeric(grepl("trt",colnames(y))))
dds <- DESeqDataSetFromMatrix(countData = y, colData = samplesheet, design = ~ trt )
res <- DESeq(dds)
z<- DESeq2::results(res)
vsd <- vst(dds, blind=FALSE)
zz<-cbind(z,assay(vsd))
x[[6]]<-as.data.frame(zz[order(zz$padj),])
sig<-subset(zz,padj<0.05)
x[[7]]<-rownames(sig[which(sig$log2FoldChange>0),])
x[[8]]<-rownames(sig[which(sig$log2FoldChange<0),])
x
}


#################################################
# define mitch function
##################################################
run_mitch<-function(y,DGE_FUNC,gsets, N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS) {
library("mitch")
dge<-sapply(y,"[",6)
names(dge)<-paste0("x",1:length(dge),sep="")
w<-mitch_import(dge, DGE_FUNC , joinType="full")

for (N in 1:ncol(w)) {
  ww<-w[,N,drop=F]
  res<-mitch_calc(ww,gsets,priority="significance",cores=8)
  y[[N]][[9]]<-res$enrichment_result[which(res$enrichment_result$s.dist>0 & res$enrichment_result$p.adjustANOVA<0.05 ),1]
  y[[N]][[10]]<-res$enrichment_result[which(res$enrichment_result$s.dist<0 & res$enrichment_result$p.adjustANOVA<0.05 ),1]
}

obs_up<-sapply(y,"[",9)
obs_dn<-sapply(y,"[",10)

gt_up<-sapply(y,"[",4)
gt_up<-lapply( gt_up , names)
gt_dn<-sapply(y,"[",5)
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
nrows<-as.numeric(lapply( sapply(y,"[",1 ), nrow))
true_neg<-mean(nrows-(true_pos+false_pos+false_neg))

p<-true_pos/(true_pos+false_pos)
r<-true_pos/(true_pos+false_neg)
f<-2*p*r/(p+r)

attr(y,'mitch_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
y
}

##################################
# hypergeometric test function (limited to deseq2)
##################################
run_hypergeometric<-function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge<-sapply(x,"[",6)

ups<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )

l_ups<-sapply(ups,length)
l_dns<-sapply(dns,length)
geneset_sizes<-sapply( gsets , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {
  universe=length(rownames(dge[[d]]))
  n_ups[[d]]<-sapply( 1:length(gsets),function(x){length(which(gsets[[x]] %in% ups[[d]] ))} )
  n_dns[[d]]<-sapply( 1:length(gsets),function(x){length(which(gsets[[x]] %in% dns[[d]] ))} )

  p_ups[[d]]<-sapply( 1:length(gsets),function(x){phyper((n_ups[[d]][[x]]-1),l_ups[[d]],universe-geneset_sizes[[x]],geneset_sizes[[x]],lower.tail=FALSE,log.p=FALSE)})
  p_dns[[d]]<-sapply( 1:length(gsets),function(x){phyper((n_dns[[d]][[x]]-1),l_dns[[d]],universe-geneset_sizes[[x]],geneset_sizes[[x]],lower.tail=FALSE,log.p=FALSE)})

  x[[d]][[11]]<-names(gsets[which(p.adjust(p_ups[[d]],method="fdr")<0.05)])
  x[[d]][[12]]<-names(gsets[which(p.adjust(p_dns[[d]],method="fdr")<0.05)])
}

obs_up<-sapply(x,"[",11)
obs_dn<-sapply(x,"[",12)

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

attr(x,'phyper_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# 1-sided fisher test function
##################################
run_fisher <- function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge<-sapply(x,"[",6)

ups<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns<-lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )
bgs <- lapply(dge, function(x) { rownames(x) })

l_ups<-sapply(ups,length)
l_dns<-sapply(dns,length)
geneset_sizes<-sapply( gsets , length )

# calculate number of genes in sets that are up and downregulated
upreg_incat=gsdet=unreg_incat=upreg_not_incat=unreg_not_incat=NULL
upreg_incat=gsdet=unreg_incat=upreg_not_incat=unreg_not_incat=list()
dnreg_incat=gsdet=unreg_incat=dnreg_not_incat=unreg_not_incat=NULL
dnreg_incat=gsdet=unreg_incat=dnreg_not_incat=unreg_not_incat=list()
p_ups=p_dns=NULL
p_ups=p_dns=list()

for (d in 1:length(dge)) {
  universe =length(rownames(dge[[d]]))

  upreg_incat[[d]]<-sapply( 1:length(gsets),function(y){length(which(gsets[[y]] %in% ups[[d]] ))} )
  gsdet[[d]] <- sapply( 1:length(gsets),function(y){length(which(gsets[[y]] %in% bgs[[d]] ))} )
  unreg_incat[[d]] <- gsdet[[d]] - upreg_incat[[d]]
  upreg_not_incat[[d]]<-sapply( 1:length(gsets),function(y){length(which(!ups[[d]] %in% gsets[[y]] ))} )
  unreg_not_incat[[d]] <- universe - upreg_incat[[d]] - unreg_incat[[d]] - upreg_not_incat[[d]]

  p_ups[[d]] <- sapply( 1:length(gsets),function(y){
    mx <- matrix(c( upreg_incat[[d]][[y]] , unreg_incat[[d]][[y]], upreg_not_incat[[d]][[y]], unreg_not_incat[[d]][[y]]),ncol=2)

    fres <- fisher.test(mx,alternative = "greater")
    fres$p
  })

  dnreg_incat[[d]]<-sapply( 1:length(gsets),function(y){length(which(gsets[[y]] %in% dns[[d]] ))} )
  gsdet[[d]] <- sapply( 1:length(gsets),function(y){length(which(gsets[[y]] %in% bgs[[d]] ))} )
  unreg_incat[[d]] <- gsdet[[d]] - dnreg_incat[[d]]
  dnreg_not_incat[[d]]<-sapply( 1:length(gsets),function(y){length(which(!dns[[d]] %in% gsets[[y]] ))} )
  unreg_not_incat[[d]] <- universe - dnreg_incat[[d]] - unreg_incat[[d]] - dnreg_not_incat[[d]]

  p_dns[[d]] <- sapply( 1:length(gsets),function(y){
    mx <- matrix(c( dnreg_incat[[d]][[y]] , unreg_incat[[d]][[y]], dnreg_not_incat[[d]][[y]], unreg_not_incat[[d]][[y]] ),ncol=2)
    fres <- fisher.test(mx,alternative = "greater")
    fres$p
  })

  x[[d]][[19]]<-names(gsets[which(p.adjust(p_ups[[d]],method="fdr")<0.05)])
  x[[d]][[20]]<-names(gsets[which(p.adjust(p_dns[[d]],method="fdr")<0.05)])
}

obs_up<-sapply(x,"[",19)
obs_dn<-sapply(x,"[",20)

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

attr(x,'fisher_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# clusterprofiler default function
##################################
# Note that clusterprofiler requires different gene set format
run_clusterprofiler_default <-function(x,DGE_FUNC,gsets2,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge <- sapply(x,"[",6)

ups <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )
bgs <- lapply(dge, function(x) { rownames(x) } )

l_ups <- sapply(ups,length)
l_dns <- sapply(dns,length)
geneset_sizes <- sapply( gsets2 , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {

# clusterprofiler UP
ora_up <- as.data.frame(enricher(gene = ups[[d]] ,
  universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_up$geneID <- NULL
ora_up <- subset(ora_up,p.adjust<0.05 )
ora_ups <- rownames(ora_up)
obs_up[[d]] <- ora_ups

# clusterprofiler DOWN
ora_dn <- as.data.frame(enricher(gene = dns[[d]] ,
  universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_dn$geneID <- NULL
ora_dn <- subset(ora_dn,p.adjust<0.05 )
ora_dns <- rownames(ora_dn)
obs_dn[[d]] <- ora_dns

x[[d]][[13]] <- ora_ups
x[[d]][[14]] <- ora_dns

}

#ground truth comparison
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

attr(x,'clusterprofiler_default') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# clusterprofiler fixed function
##################################
# Note that clusterprofiler requires different gene set format
run_clusterprofiler_fixed <-function(x,DGE_FUNC,gsets2,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge <- sapply(x,"[",6)

ups <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )
bgs <- lapply(dge, function(x) { rownames(x) } )

l_ups <- sapply(ups,length)
l_dns <- sapply(dns,length)
geneset_sizes <- sapply( gsets2 , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {

bgset <- bgs[[d]]
bgdf <- data.frame(term="bglist",gene=bgset)
gsets3 <- rbind(gsets2,bgdf)

# clusterprofiler UP
ora_up <- as.data.frame(enricher(gene = ups[[d]] ,
  universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets3,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1 ))

ora_up$geneID <- NULL
ora_up <- subset(ora_up,p.adjust<0.05 & Count)
ora_ups <- rownames(ora_up)
obs_up[[d]] <- ora_ups

# clusterprofiler DOWN
ora_dn <- as.data.frame(enricher(gene = dns[[d]] ,
  universe = bgs[[d]],  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets3,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1 ))

ora_dn$geneID <- NULL
ora_dn <- subset(ora_dn,p.adjust<0.05 & Count)
ora_dns <- rownames(ora_dn)
obs_dn[[d]] <- ora_dns

x[[d]][[15]] <- ora_ups
x[[d]][[16]] <- ora_dns

}

#ground truth comparison
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

attr(x,'clusterprofiler_fixed') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}


##################################
# clusterprofiler nobg
##################################
run_clusterprofiler_nobg <-function(x,DGE_FUNC,gsets2,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge <- sapply(x,"[",6)

ups <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )
bgs <- lapply(dge, function(x) { rownames(x) } )

l_ups <- sapply(ups,length)
l_dns <- sapply(dns,length)
geneset_sizes <- sapply( gsets2 , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {

# clusterprofiler UP
ora_up <- as.data.frame(enricher(gene = ups[[d]] ,
  universe = NULL,  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_up$geneID <- NULL
ora_up <- subset(ora_up,p.adjust<0.05 )
ora_ups <- rownames(ora_up)
obs_up[[d]] <- ora_ups

# clusterprofiler DOWN
ora_dn <- as.data.frame(enricher(gene = dns[[d]] ,
  universe = NULL,  minGSSize=2, maxGSSize = 500000, TERM2GENE = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_dn$geneID <- NULL
ora_dn <- subset(ora_dn,p.adjust<0.05 )
ora_dns <- rownames(ora_dn)
obs_dn[[d]] <- ora_dns

x[[d]][[17]] <- ora_ups
x[[d]][[18]] <- ora_dns

}

#ground truth comparison
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

attr(x,'clusterprofiler_nobg') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# DOSE enrich internal
##################################
EXTID2TERMID <- function(gene, USER_DATA) {
    if (inherits(USER_DATA, "environment")) { 
        EXTID2PATHID <- get("EXTID2PATHID", envir = USER_DATA)

        qExtID2Path <- EXTID2PATHID[gene]
    } else if (inherits(USER_DATA, "GSON")) {
        gsid2gene <- USER_DATA@gsid2gene
        qExtID2Path <- setNames(lapply(gene, function(x) {
            subset(gsid2gene, gsid2gene$gene == x)[["gsid"]]
        }), gene)
    } else {
        stop("not supported")
    }

    len <- sapply(qExtID2Path, length)
    notZero.idx <- len != 0
    qExtID2Path <- qExtID2Path[notZero.idx]

    return(qExtID2Path)
}

ALLEXTID <- function(USER_DATA) {
    if (inherits(USER_DATA, "environment")) { 
        PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
        res <- unique(unlist(PATHID2EXTID))
    } else if (inherits(USER_DATA, "GSON")) {
        gsid2gene <- USER_DATA@gsid2gene
        res <- unique(gsid2gene$gene)
    } else {
        stop("not supported")
    }

    return(res)
}

TERMID2EXTID <- function(term, USER_DATA) {
    if (inherits(USER_DATA, "environment")) { 
        PATHID2EXTID <- get("PATHID2EXTID", envir = USER_DATA)
        res <- PATHID2EXTID[term]
    } else if (inherits(USER_DATA, "GSON")) {
        gsid2gene <- USER_DATA@gsid2gene
        res <- setNames(lapply(term, function(x) {
            subset(gsid2gene, gsid2gene$gsid == x)[["gene"]]
        }), term)
    } else {
        stop("not supported")
    }    

    return(res)
}

get_geneSet_index <- function(geneSets, minGSSize, maxGSSize) {
    if (is.na(minGSSize) || is.null(minGSSize))
        minGSSize <- 1
    if (is.na(maxGSSize) || is.null(maxGSSize))
        maxGSSize <- Inf #.Machine$integer.max

    ## index of geneSets in used.
    ## logical
    geneSet_size <- sapply(geneSets, length)
    idx <-  minGSSize <= geneSet_size & geneSet_size <= maxGSSize
    return(idx)
}

TERM2NAME <- function(term, USER_DATA) {
    if (inherits(USER_DATA, "environment")) { 
        PATHID2NAME <- get("PATHID2NAME", envir = USER_DATA)
        #if (is.null(PATHID2NAME) || is.na(PATHID2NAME)) {
        if (is.null(PATHID2NAME) || all(is.na(PATHID2NAME))) {
            return(as.character(term))
        }
        return(PATHID2NAME[term])
    } else if (inherits(USER_DATA, "GSON")) {
        gsid2name <- USER_DATA@gsid2name
        res <- setNames(vapply(term, function(x) {
            subset(gsid2name, gsid2name$gsid == x)[["name"]]
        }, character(1)), term)
        return(res)
    } 

    return(as.character(term)) 
}

enricher_internal <- function(gene,
                              pvalueCutoff,
                              pAdjustMethod="BH",
                              universe = NULL,
                              minGSSize=10,
                              maxGSSize=500,
                              qvalueCutoff=0.2,
                              USER_DATA){

    ## query external ID to Term ID
    gene <- as.character(unique(gene))
    qExtID2TermID <- EXTID2TERMID(gene, USER_DATA)
    qTermID <- unlist(qExtID2TermID)
    if (is.null(qTermID)) {
        message("--> No gene can be mapped....")
        if (inherits(USER_DATA, "environment")) {
            p2e <- get("PATHID2EXTID", envir=USER_DATA)
            sg <- unique(unlist(p2e[1:10]))
        } else {
            sg <- unique(USER_DATA@gsid2gene$gene[1:100])
        }
        sg <- sample(sg, min(length(sg), 6))
        message("--> Expected input gene ID: ", paste0(sg, collapse=','))

        message("--> return NULL...")
        return(NULL)
    }

    ## Term ID -- query external ID association list.
    qExtID2TermID.df <- data.frame(extID=rep(names(qExtID2TermID),
                                             times=lapply(qExtID2TermID, length)),
                                   termID=qTermID)
    qExtID2TermID.df <- unique(qExtID2TermID.df)

    qTermID2ExtID <- with(qExtID2TermID.df,
                          split(as.character(extID), as.character(termID)))

    extID <- ALLEXTID(USER_DATA)
    if (missing(universe))
        universe <- NULL
    if(!is.null(universe)) {
        if (is.character(universe)) {
            force_universe <- getOption("enrichment_force_universe", FALSE)
            if (!force_universe) {
                extID <- intersect(extID, universe)
            }
        } else {
            ## https://github.com/YuLab-SMU/clusterProfiler/issues/217
            message("`universe` is not in character and will be ignored...")
        }
    }

    qTermID2ExtID <- lapply(qTermID2ExtID, intersect, extID)

    ## Term ID annotate query external ID
    qTermID <- unique(names(qTermID2ExtID))


    termID2ExtID <- TERMID2EXTID(qTermID, USER_DATA)
    termID2ExtID <- lapply(termID2ExtID, intersect, extID)

    geneSets <- termID2ExtID

    idx <- get_geneSet_index(termID2ExtID, minGSSize, maxGSSize)

    if (sum(idx) == 0) {
        msg <- paste("No gene sets have size between", minGSSize, "and", maxGSSize, "...")
        message(msg)
        message("--> return NULL...")
        return (NULL)
    }

    termID2ExtID <- termID2ExtID[idx]
    qTermID2ExtID <- qTermID2ExtID[idx]
    qTermID <- unique(names(qTermID2ExtID))

    ## prepare parameter for hypergeometric test
    k <- sapply(qTermID2ExtID, length)
    k <- k[qTermID]
    M <- sapply(termID2ExtID, length)
    M <- M[qTermID]

    N <- rep(length(extID), length(M))
    ## n <- rep(length(gene), length(M)) ## those genes that have no annotation should drop.
    n <- rep(length(qExtID2TermID), length(M))
    args.df <- data.frame(numWdrawn=k-1, ## White balls drawn
                          numW=M,        ## White balls
                          numB=N-M,      ## Black balls
                          numDrawn=n)    ## balls drawn


    ## calcute pvalues based on hypergeometric model
    pvalues <- apply(args.df, 1, function(n)
                     phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
                     )

    ## gene ratio and background ratio
    GeneRatio <- apply(data.frame(a=k, b=n), 1, function(x)
                       paste(x[1], "/", x[2], sep="", collapse="")
                       )
    BgRatio <- apply(data.frame(a=M, b=N), 1, function(x)
                     paste(x[1], "/", x[2], sep="", collapse="")
                     )


    Over <- data.frame(ID = as.character(qTermID),
                       GeneRatio = GeneRatio,
                       BgRatio = BgRatio,
                       pvalue = pvalues,
                       stringsAsFactors = FALSE)

    p.adj <- p.adjust(Over$pvalue, method=pAdjustMethod)
    qobj <- tryCatch(qvalue(p=Over$pvalue, lambda=0.05, pi0.method="bootstrap"), error=function(e) NULL)

    # if (class(qobj) == "qvalue") {
    if (inherits(qobj, "qvalue")) {
        qvalues <- qobj$qvalues
    } else {
        qvalues <- NA
    }

    geneID <- sapply(qTermID2ExtID, function(i) paste(i, collapse="/"))
    geneID <- geneID[qTermID]
    Over <- data.frame(Over,
                       p.adjust = p.adj,
                       qvalue = qvalues,
                       geneID = geneID,
                       Count = k,
                       stringsAsFactors = FALSE)

    Description <- TERM2NAME(qTermID, USER_DATA)

    if (length(qTermID) != length(Description)) {
        idx <- qTermID %in% names(Description)
        Over <- Over[idx,]
    }
    Over$Description <- Description
    nc <- ncol(Over)
    Over <- Over[, c(1,nc, 2:(nc-1))]


    Over <- Over[order(pvalues),]


    Over$ID <- as.character(Over$ID)
    Over$Description <- as.character(Over$Description)

    row.names(Over) <- as.character(Over$ID)

    x <- new("enrichResult",
             result         = Over,
             pvalueCutoff   = pvalueCutoff,
             pAdjustMethod  = pAdjustMethod,
             qvalueCutoff   = qvalueCutoff,
             gene           = as.character(gene),
             universe       = extID,
             geneSets       = geneSets,
             organism       = "UNKNOWN",
             keytype        = "UNKNOWN",
             ontology       = "UNKNOWN",
             readable       = FALSE
             )
    if (inherits(USER_DATA, "GSON")) {
        if (!is.null(USER_DATA@keytype)) {
            x@keytype <- USER_DATA@keytype
        }
        if (!is.null(USER_DATA@species)) {
            x@organism <- USER_DATA@species
        } 
        if (!is.null(USER_DATA@gsname)) {
            x@ontology <- gsub(".*;", "", USER_DATA@gsname)
        }
    }
    return (x)
}

run_dose <-function(x,DGE_FUNC,gsets2,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge <- sapply(x,"[",6)

ups <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange > 0)) } )
dns <- lapply(dge, function(x) { rownames(subset(x, padj<0.05 & log2FoldChange < 0)) } )
bgs <- lapply(dge, function(x) { rownames(x) } )

l_ups <- sapply(ups,length)
l_dns <- sapply(dns,length)
geneset_sizes <- sapply( gsets2 , length )

# calculate number of genes in sets that are up and downregulated
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=NULL
n_dns=n_ups=p_ups=p_dns=obs_up=obs_dn=list()

for (d in 1:length(dge)) {

# DOSE UP
ora_up <- as.data.frame(enricher_internal(gene = ups[[d]] ,
  universe = NULL,  minGSSize=2, maxGSSize = 500000, USER_DATA = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_up$geneID <- NULL
ora_up <- subset(ora_up,p.adjust<0.05 )
ora_ups <- rownames(ora_up)
obs_up[[d]] <- ora_ups

# DOSE DOWN
ora_dn <- as.data.frame(enricher_internal(gene = dns[[d]] ,
  universe = NULL,  minGSSize=2, maxGSSize = 500000, USER_DATA = gsets2,
  pAdjustMethod="fdr",  pvalueCutoff = 1, qvalueCutoff = 1  ))

ora_dn$geneID <- NULL
ora_dn <- subset(ora_dn,p.adjust<0.05 )
ora_dns <- rownames(ora_dn)
obs_dn[[d]] <- ora_dns

x[[d]][[21]] <- ora_ups
x[[d]][[22]] <- ora_dns

}

#ground truth comparison
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

attr(x,'dose_enricher') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x

}

##################################
# FGSEA function
##################################
run_fgsea<-function(x,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS){

dge<-sapply(x,"[",6)

xx<-lapply( dge , function(x) { 
 s<-x$stat
 names(s)<-rownames(x)
 p<-as.data.frame(fgsea(pathways=gsets, stats=s ))
 p
} )

obs_up<-lapply(xx, function(x) { subset(x,padj<0.05 & ES>0)[,1] } )
obs_dn<-lapply(xx, function(x) { subset(x,padj<0.05 & ES<0)[,1] } )

for (d in 1:length(dge)) {
  x[[d]][[9]]<-obs_up[[d]]
  x[[d]][[10]]<-obs_dn[[d]]
}

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

attr(x,'fgsea_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
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

attr(x,'gst_res') <-data.frame(N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,true_pos,false_pos,true_neg,false_neg,p,r,f)
x
}


##################################
# aggregate function
##################################
agg_dge<-function(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets) {

#TEST# N_REPS=5 ; SUM_COUNT=30000000 ; VARIANCE=0.2 ; FRAC_DE=0.05 ; FC=1 ; SIMS=8 ; DGE_FUNC="deseq2" ; gsets=gsets

xxx <- RepParallel(SIMS,simrna(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,gsets), simplify=F, mc.cores = 8 )

# run deseq2
xxx <- mclapply(xxx , DGE_FUNC , mc.cores = 8 )

# run clusterprofiler default
writeGmtPathways(pathways=gsets,gmt.file="mypathways.gmt")
gsets2 <- read.gmt("mypathways.gmt")
xxx <- run_clusterprofiler_default(x=xxx,DGE_FUNC,gsets2,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run clusterprofiler fixed
xxx <- run_clusterprofiler_fixed(x=xxx,DGE_FUNC,gsets2,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run clusterprofiler nobg
xxx <- run_clusterprofiler_nobg(x=xxx,DGE_FUNC,gsets2,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run dose
#xxx <- run_dose(x=xxx,DGE_FUNC,gsets2,N_REPS=N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run phyper
xxx<-run_hypergeometric(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run fgsea
xxx<-run_fgsea(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

# run fisher
xxx<-run_fisher(xxx,DGE_FUNC,gsets,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS)

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

