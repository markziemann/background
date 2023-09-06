library("tidyverse")
library("parallel")
library("edgeR")
library("DESeq2")
library("limma")
library("stringi")

source("simpw_func.R")

# obtain count data
a<-countData()

# generate some gene sets
gsets<-randomGeneSets(a)


# test
# x<-agg_dge(a,3,40000000,0.3,0.05,1,10,"deseq2",gsets) ; str(x)

###############################################
# run simulations over a range of parameters
###############################################
SIMS=10
unlink("simpw_res_running.tsv")
res=NULL
for ( FRAC_DE in c(0.05)) {
  for (FC in c(1)) {
    for (N_REPS in c(5)) {
      for (DGE_FUNC in c("deseq2")) {
        for ( SUM_COUNT in c(50000000)) {
          for ( VARIANCE in c(0.5)) {
            x<-agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
            x<-as.data.frame(do.call(rbind, x))
            write.table(x,file="simpw_res_running.tsv",quote=F,sep='\t',append=T)
            res=c(res,x)
          }
        }
      }
    }
  }
}
