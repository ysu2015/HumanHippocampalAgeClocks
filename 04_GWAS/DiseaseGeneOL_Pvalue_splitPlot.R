#library('Seurat')
packageVersion("Seurat") # 4.0.5
library(dplyr)

GWASDir <- "/media/yijingsu/BackUp03/00_2023_Neuron_SPLiTSeq_HumanHippocampus_MolecularClocks_2022_0701_To_2023_1030/01_ScRNAseq_Analysis/01_HIPP_ClockGenes/04_GWAS/"
################################
setwd("/media/yijingsu/BackUp03/00_2023_Neuron_SPLiTSeq_HumanHippocampus_MolecularClocks_2022_0701_To_2023_1030/02_afterFirstSubmission/31_PosNegClockGenesGWAS/")

risk.ASD<-read.table(paste(GWASDir,"./MAGMA_10Up15Down/ASD_iPSYCH-PGC_ASD_Nov2017.10UP.1.5DOWN.genes.out", sep=""), header = T)
head(risk.ASD)

risk.ASD.sig<-filter(risk.ASD, P < 0.05)

risk.ad = read.table(paste(GWASDir,"./MAGMA_10Up15Down/AD.10UP.1.5DOWN.genes.out",  sep=""),header = T)
risk.ad.sig = filter(risk.ad, P < 0.05)

risk.epi = read.table(paste(GWASDir,"./MAGMA_10Up15Down/all_epilepsy_METAL.10UP.1.5DOWN.genes.out",  sep=""),header = T)
risk.epi.sig = filter(risk.epi, P < 0.05)

risk.scz = read.table(paste(GWASDir,"./MAGMA_10Up15Down/PGC3_SCZ_wave3.european.autosome.public.v3_RemovedFirst74.10UP.1.5DOWN.genes.out", sep=""), header = T)
risk.scz.sig = filter(risk.scz, P < 0.05)

risk.Bip<-read.table(paste(GWASDir,"./MAGMA_10Up15Down/BP_pgc-bip2021-all_RemovedFirst73_LastColumn.10UP.1.5DOWN.genes.out", sep=""), header = T)
risk.Bip.sig<-filter(risk.Bip, P < 0.05)

risk.Anxiety<-read.table(paste(GWASDir,"./MAGMA_10Up15Down/anxiety.meta.full.cc.10UP.1.5DOWN.genes.out",  sep=""),header = T)
risk.Anxiety.sig<-filter(risk.Anxiety, P < 0.05)

risk.Depressive<-read.table(paste(GWASDir,"./MAGMA_10Up15Down/MDD_jamapsy_Giannakopoulou_2021.10UP.1.5DOWN.genes.out",  sep=""),header = T)
risk.Depressive.sig<-filter(risk.Depressive, P < 0.05)

###
Disease<-list("risk.ASD"=risk.ASD.sig,
              "risk.ad"=risk.ad.sig,
              "risk.epi"=risk.epi.sig,
              "risk.scz"=risk.scz.sig,
              "risk.Bip"=risk.Bip.sig,
              "risk.Anxiety"=risk.Anxiety.sig,
              #   "risk.Dementia"=risk.Dementia,
              "risk.Depressive"=risk.Depressive.sig
              #   "risk.DevDisabilities"=risk.DevDisabilities,
              #    "risk.Mood"=risk.Mood,
              #   "risk.Retardation"=risk.Retardation,
              #    "Mental.Disorders"=Mental.Disorders,
              #    "PD"=PD
)
Diseasename<-names(Disease)
####################
#########################
Total.clock.gene<-read.csv("ClockGeneCutoff.csv")
head(Total.clock.gene)
dim(Total.clock.gene) #3296


# not cell type specific, include shared genes ##
gene.list.OPC<-filter(Total.clock.gene, Cell.types =="OPC");
  dim(gene.list.OPC) #990
  pos.gene.list.OPC <- filter(gene.list.OPC, Correlations == "positive"); dim(pos.gene.list.OPC) #547
  neg.gene.list.OPC <- filter(gene.list.OPC, Correlations == "negative"); dim(neg.gene.list.OPC) #443
  #pos.gene.list.OPC <- filter(gene.list.OPC, Coefficients > 0); dim(pos.gene.list.OPC) #547
  #neg.gene.list.OPC <- filter(gene.list.OPC, Coefficients < 0 ); dim(neg.gene.list.OPC) #443
gene.list.OLIG<-filter(Total.clock.gene, Cell.types =="OLIG");
  dim(gene.list.OLIG) #597
  pos.gene.list.OLIG <- filter(gene.list.OLIG, Correlations == "positive"); dim(pos.gene.list.OLIG) #358
  neg.gene.list.OLIG <- filter(gene.list.OLIG, Correlations == "negative"); dim(neg.gene.list.OLIG) #239
 #  pos.gene.list.OLIG <- filter(gene.list.OLIG, Coefficients > 0); dim(pos.gene.list.OLIG) #547
  # neg.gene.list.OLIG <- filter(gene.list.OLIG, Coefficients < 0 ); dim(neg.gene.list.OLIG) #443
gene.list.AST<-filter(Total.clock.gene, Cell.types =="AST");
   dim(gene.list.AST) #679
   pos.gene.list.AST <- filter(gene.list.AST, Correlations == "positive"); dim(pos.gene.list.AST) #371
   neg.gene.list.AST <- filter(gene.list.AST, Correlations == "negative"); dim(neg.gene.list.AST) #308
  # pos.gene.list.AST <- filter(gene.list.AST, Coefficients > 0); dim(pos.gene.list.AST) #371
  # neg.gene.list.AST <- filter(gene.list.AST, Coefficients < 0 ); dim(neg.gene.list.AST) #308
gene.list.AST<-filter(Total.clock.gene, Cell.types =="AST");
   dim(gene.list.AST) #679
   pos.gene.list.AST <- filter(gene.list.AST, Correlations == "positive"); dim(pos.gene.list.AST) #371
   neg.gene.list.AST <- filter(gene.list.AST, Correlations == "negative"); dim(neg.gene.list.AST) #308
   # pos.gene.list.AST <- filter(gene.list.AST, Coefficients > 0); dim(pos.gene.list.AST) #371
   # neg.gene.list.AST <- filter(gene.list.AST, Coefficients < 0 ); dim(neg.gene.list.AST) #308   
gene.list.MG<-filter(Total.clock.gene, Cell.types =="MG");
   dim(gene.list.MG) #1030
   pos.gene.list.MG <- filter(gene.list.MG, Correlations == "positive"); dim(pos.gene.list.MG) #595
   neg.gene.list.MG <- filter(gene.list.MG, Correlations == "negative"); dim(neg.gene.list.MG) #435
    #pos.gene.list.MG <- filter(gene.list.MG, Coefficients > 0); dim(pos.gene.list.MG) #371
    #neg.gene.list.MG <- filter(gene.list.MG, Coefficients < 0 ); dim(neg.gene.list.MG) #308 
library("GeneOverlap")
library("org.Hs.eg.db")
library(annotate)
library(dplyr) 
   
#CandidateGene <- pos.gene.list.MG$Gene
CandidateGene <- neg.gene.list.AST$Gene

results_list <- list()     
for(i in Diseasename)
  {print(i)
   
    DiseaseGene<-Disease[i][[1]]$GENE
    DiseaseGene.SYMBOL = getSYMBOL(as.character(DiseaseGene),data='org.Hs.eg')
    go.obj<-newGeneOverlap(CandidateGene,DiseaseGene.SYMBOL,genome.size=23936) #
  #   go.obj
    go.obj <- testGeneOverlap(go.obj)
   # print(go.obj)
    overlap_result <- data.frame(
      # jaccard_index = getJaccardIndex(go.obj),
      p_value = getPval(go.obj),
      odds_ratio = getOddsRatio(go.obj)
    )
    # Store each result
    results_list[[i]] <- overlap_result
}
  
# Combine all results into one data.frame
final_results <- do.call(rbind, results_list)   
   
# View results
print(final_results)

''' output different cell types p value to I:/32Samples/DiseaseAssociation/DiseasePvalue.csv for plot '''

  ##########################  
library(dplyr) 
  all.GO<-read.csv("DiseasePvalue.csv")
  head(all.GO)
  dim(all.GO)
  summary(all.GO)
  
  colnames(all.GO)<-c("Correlation","CellType","DiseaseCat","pValue","Odds_ratio")
  
 # GO<-filter(all.GO, Correlation == "Pos")
  GO <- filter(all.GO, Correlation == "Neg")
  dim(GO)
  

  library(ggplot2)
  #jpeg("ASD.Disease.jpeg",width = 800,height = 400)
  
ggplot(GO, aes(x=factor(CellType,
                               level=c("AST","OPC","OLIG","MG")),
                      y=factor(DiseaseCat,
                               level=rev(c("ASD","AD","SCZ","EPI","BP","Anxiety","MDD")))))+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_point(aes(size=Odds_ratio, col=-log10(pValue)))+
    scale_size_continuous(range = c(1,5)) +
    scale_color_continuous(low = "grey", high = "red2")

ggsave("PvalueDisease_pos.pdf",width=3.2,height=3)
ggsave("PvalueDisease_neg.pdf",width=3.2,height=3)
range(-log10(GO$pValue))
range(GO$pValue)
  #dev.off()
  ##########################  ##########################  ##########################  ##########################    
##########################  ##########################  ##########################  ##########################  
##########################  ##########################  ##########################  ##########################  
##########################  ##########################  ##########################  ##########################   
