# Convert GMT to RDS
library("clusterProfiler")
l <- list.files(".",pattern = "gmt")
l <- l[order(l)]
gs <- lapply(l,function(x) {
  y <- read.gmt(x)
  y$term <- gsub("_"," ",y$term)
  return(y)
})
str(gs)
names(gs) <- c("KEGG Pathways",
               "Reactome Pathways",
               "WikiPathways",
               "miR targets",
               "GTRD TF Targets", 
               "Gene Ontology",
               "Human Phenotype Ontology",
               "CellMarkers", 
               "Hallmark")
                            
saveRDS(gs,"gs.Rds")
