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
gsets<-randomGeneSets(a,setsize=30,nsets=500)

# test
# x<-agg_dge(a=a,N_REPS=5,SUM_COUNT=50000000,VARIANCE=0.4,FRAC_DE=0.1,FC=1,SIMS=8,"deseq2",gsets)
# str(x)

###############################################
# run simulations over a range of parameters
###############################################
SIMS=8
FRAC_DE=0.05
FC=1
N_REPS=5
DGE_FUNC="deseq2"
SUM_COUNT=3e7
VARIANCE=c(0.2,0.3,0.4,0.5)

mygrid <- expand.grid(FRAC_DE,FC,N_REPS,DGE_FUNC,SUM_COUNT,VARIANCE)
colnames(mygrid) <- c("FRAC_DE","FC","N_REPS","DGE_FUNC","SUM_COUNT","VARIANCE")

mygrid

res <- lapply(1:nrow(mygrid), function(i) {
  FRAC_DE=mygrid[i,"FRAC_DE"]
  FC=mygrid[i,"FC"]
  N_REPS=mygrid[i,"N_REPS"]
  DGE_FUNC=as.character(mygrid[i,"DGE_FUNC"])
  SUM_COUNT=mygrid[i,"SUM_COUNT"]
  VARIANCE=mygrid[i,"VARIANCE"]
  x <- agg_dge(a,N_REPS,SUM_COUNT,VARIANCE,FRAC_DE,FC,SIMS,DGE_FUNC,gsets)
  as.data.frame(do.call(rbind, x))
})

res
