# Train with full set of data. Tidy form.

library(tidyverse)
library(glmnet)
library(caret)

setwd("./")

CT<- c("AST","OLIG","OPC","MG")

#CT<- "OPC"
for (Celltypes in CT) {
  print(Celltypes)
df <- readRDS(paste0("Models/Model_bootstrap_pseudocell",Celltypes,"_15_seed42.rds",sep=""))
#df <- df %>% ungroup %>% select(-c(hash.ID, orig.ident, data))
df <- df %>% ungroup %>% select(-c(hash.ID, orig.ident))
df
lognorm <- function(input) {
    norm <- sweep(input, MARGIN = 1, FUN = "/", STATS = rowSums(input))
    log1p(norm * 10000)
}

df2 <- df %>% mutate(lognorm = map(pseudocell_all, lognorm))
df2$pseudocell_all <- NULL
df2 <- unnest(df2, lognorm)
df <- NULL # free up memory

by_celltype <- df2 %>% group_by(Celltype.LowRes) %>% nest()
colnames(by_celltype)[2] <- "lognormalized"

# Model function
celltype_model <- function(input) {
    cv.glmnet(x = as.matrix(input[, -1]), y = as.matrix(input[, 1]),
              type.measure = "mae", standardize = F, relax = F, nfolds = 5)
}

# Apply model function using map()
models <- by_celltype %>% mutate(model = map(lognormalized, celltype_model))

models
#   Celltype.LowRes lognormalized             model     
#   <fct>           <list>                    <list>    
# 1 Oligodendro     <tibble [2,800 x 20,949]> <cv.glmnt>
# 2 Endothelial     <tibble [2,800 x 20,949]> <cv.glmnt>
# 3 Microglia       <tibble [2,800 x 20,949]> <cv.glmnt>
# 4 Astrocyte_qNSC  <tibble [2,800 x 20,949]> <cv.glmnt>
# 5 aNSC_NPC        <tibble [2,800 x 20,949]> <cv.glmnt>
# 6 Neuroblast      <tibble [2,800 x 20,949]> <cv.glmnt>

saveRDS(models, paste0("./Models/models_",Celltypes,"_bootstrap.rds",sep=""))
}
