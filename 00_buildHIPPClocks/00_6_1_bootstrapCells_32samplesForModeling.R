
rm(list=ls())
library(Seurat)
library(tidyverse)
library(ape)

setwd("./")

all<-readRDS("32sample.4GliaTypes.Down1500.SCTV1.counts.rds")
#===================================================================================================


CT<- c("AST","OLIG","Endothelial","OPC","MG","GABA.N","Ependymal")#, "Cajal-Retzius","Choroid.plexus")#"Glut.N",

#==================================================================================================
#==================================================================================================
#==================================================================================================

## normalization of RNA counts within cell level
## svz1<-NormalizeData(svz1, normalization.method = "LogNormalize", scale.factor = 10000)

for (Celltypes in CT) {
  print(Celltypes)
svz <- subset(all, subset = Celltype.LowRes %in% Celltypes)
svz
#length(svz@meta.data$Age.ac.Year)#19589

## select for genes specific going up/down along age within each cell types
## so far using the lmfit, will switch to Moran'I later

gene.name<-read.csv(paste0("./Models/Gene.List.",Celltypes,".csv",sep=""))

################################################

Convert_to_Dataframe <- function(svz) {
    DefaultAssay(svz) <- "RNA"
   # svz[["SCT"]] <- NULL
   # svz[["LMO"]] <- NULL
   #  svz[["integrated"]]<-NULL
     Celltypes<-Celltypes
     svz <- subset(svz, subset = Celltype.LowRes %in% Celltypes)
    meta <- svz@meta.data
    meta <- meta[, c("hash.ID", "Age", "Celltype.LowRes", "orig.ident")]
    raw_counts <- t(as.matrix(svz[["RNA"]]@counts))
    raw_counts <- raw_counts[, colSums(raw_counts) > 0]
    raw_counts <- raw_counts[, colnames(raw_counts) %in% gene.name$V1] # select for lm genes
    df <- as_tibble(cbind(meta, raw_counts))
    return(df)
}

df <- Convert_to_Dataframe(svz) %>%
        group_by(Celltype.LowRes, Age, orig.ident, hash.ID) %>%
        nest()
head(df) # note: Age need to be "dbl" for further analysis
# A tibble: 6 × 5
# Groups:   Celltype.LowRes, Age, orig.ident, hash.ID [6]
#hash.ID   Age Celltype.LowRes orig.ident data                    
#<chr>   <dbl> <chr>           <chr>      <list>                  
#  1 S3       1.16 OLIG            B15        <tibble [1,286 × 1,763]>
#  2 S8       6    OLIG            B15        <tibble [525 × 1,763]>  
#  3 S4       2.08 OLIG            B15        <tibble [636 × 1,763]>  
 # 4 S19     88    OLIG            B15        <tibble [1,851 × 1,763]>
#  5 S44      0.1  OLIG            BB         <tibble [2 × 1,763]>    
 # 6 S39      0.15 OLIG            B16        <tibble [649 × 1,763]>  


#===================================================================================================
# test different number for pseudocell as
start=15
end=15
interval=5

for (i in seq(start, end, interval))
{
Size=i # set as 15 original
bootstrap.pseudocells <- function(df, size=Size, n=100, replace="dynamic") {
  pseudocells <- c()
  # If dynamic then only sample with replacement if required due to shortage of cells.
  if (replace == "dynamic") {
    if (nrow(df) <= size) {replace <- TRUE} else {replace <- FALSE}
  }
  for (i in c(1:n)) {
    batch <- df[sample(1:nrow(df), size = size, replace = replace), ]
    pseudocells <- rbind(pseudocells, colSums(batch))
  }
  colnames(pseudocells) <- colnames(df)
  return(as_tibble(pseudocells))
}
    
# Apply boostrap.pseudocells using map()
set.seed(42)
df2 <- df %>% mutate(pseudocell_all = map(data, bootstrap.pseudocells)) # ~15 minutes on a laptop

#==================================================================================================
# Remove single cell data; keep just key metadata and pseudocells
df2$data <- NULL
saveRDS(df2, paste0("Models/Model_bootstrap_pseudocell",Celltypes,"_",Size,"_seed42.rds",sep=""))
print(Size)
}
}
