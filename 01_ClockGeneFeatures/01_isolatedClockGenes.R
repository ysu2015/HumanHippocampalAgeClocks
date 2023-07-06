
library(tidyverse)
library(UpSetR)
library(glmnet)
library(viridis)
library(ComplexHeatmap)

setwd("E:/Yijing/AgeClocks/12_HIPP_ClockGenes//01_isolatedClockGenes_pvalue//")

celltypes = c("OPC", "AST","MG","OLIG")

#======================================================================================================================
# Save all clock genes
all <- c()
for (celltype in celltypes) {
    print(celltype)
    models <- readRDS(paste0("./Models/models_",celltype,"_bootstrap.rds",sep=""))
    models$lognormalized <- NULL
    clock <- filter(models, Celltype.LowRes == celltype)[1,2][[1]][[1]]
    clock
        # Coefficients
    betas <- coef(clock, s = clock$lambda.min) # includes intercept, length 13194
    betas
    genes <- unlist(betas@Dimnames[1]) # includes intercept, length 13194
  #  clock.genes <- genes[which(betas !=0)]
    clock.genes <- genes[which(abs(betas) > 0.15)]
    print(paste0("Genes in model: ", length(clock.genes)))
    clock.betas <- betas[which(abs(betas) > 0.15)]
    clock.df <- data.frame("gene"=clock.genes[-1], "coef"=clock.betas[-1], "Celltype" = celltype)
    clock.df <- clock.df[order(clock.df$coef, decreasing=T),]
    write.csv(clock.df, file = paste0("./Output/clock_genes_cutoff_", celltype, ".csv"), quote = F)
    all <- rbind(all, clock.df)
}
saveRDS(all, "./Output/clock_genes_all_cutoff.rds")
