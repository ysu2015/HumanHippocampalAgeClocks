library('Seurat')
packageVersion("Seurat") # 4.0.5
library(dplyr)

################################
setwd("E:/Yijing/AgeClocks/12_HIPP_ClockGenes/04_GWAS/")

  risk.ASD<-read.table("./MAGMA_10Up15Down/ASD_iPSYCH-PGC_ASD_Nov2017.10UP.1.5DOWN.genes.out", header = T)
  head(risk.ASD)
risk.ASD.sig<-filter(risk.ASD, P < 0.05)

risk.ad = read.table("./MAGMA_10Up15Down/AD.10UP.1.5DOWN.genes.out", header = T)
risk.ad.sig = filter(risk.ad, P < 0.05)

risk.epi = read.table("./MAGMA_10Up15Down/all_epilepsy_METAL.10UP.1.5DOWN.genes.out", header = T)
risk.epi.sig = filter(risk.epi, P < 0.05)

risk.scz = read.table("./MAGMA_10Up15Down/PGC3_SCZ_wave3.european.autosome.public.v3_RemovedFirst74.10UP.1.5DOWN.genes.out", header = T)
risk.scz.sig = filter(risk.scz, P < 0.05)

risk.Bip<-read.table("./MAGMA_10Up15Down/BP_pgc-bip2021-all_RemovedFirst73_LastColumn.10UP.1.5DOWN.genes.out", header = T)
risk.Bip.sig<-filter(risk.Bip, P < 0.05)

risk.Anxiety<-read.table("./MAGMA_10Up15Down/anxiety.meta.full.cc.10UP.1.5DOWN.genes.out", header = T)
risk.Anxiety.sig<-filter(risk.Anxiety, P < 0.05)

risk.Depressive<-read.table("./MAGMA_10Up15Down/MDD_jamapsy_Giannakopoulou_2021.10UP.1.5DOWN.genes.out", header = T)
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
Total.clock.gene<-read.csv("clockGene_Table_cutoff.csv",row.names = 1)
head(Total.clock.gene)
Total.clock.gene$gene<-row.names(Total.clock.gene)




# not cell type specific, include shared genes ##
gene.list.OPC<-filter(Total.clock.gene, OPC ==1 ); gene.list.OPC$celltypes<-"OPC"
gene.list.OLIG<-filter(Total.clock.gene,  OLIG==1 ); gene.list.OLIG$celltypes<-"OLIG"
gene.list.AST<-filter(Total.clock.gene, AST==1 ); gene.list.AST$celltypes<-"AST"
gene.list.MG<-filter(Total.clock.gene, MG==1); gene.list.MG$celltypes<-"MG"

library("GeneOverlap")
library("org.Hs.eg.db")
library(annotate)

for(i in Diseasename)
  {print(i)
   
    DiseaseGene<-Disease[i][[1]]$GENE
    
    DiseaseGene.SYMBOL = getSYMBOL(as.character(DiseaseGene),data='org.Hs.eg')
   
    go.obj<-newGeneOverlap(gene.list.MG$gene,DiseaseGene.SYMBOL,genome.size=23936) #
 #   go.obj
    go.obj <- testGeneOverlap(go.obj)
    print(go.obj)
    print("\n")
    print("\n")
    print("\n")
    print("\n")
  }
''' output different cell types p value to I:/32Samples/DiseaseAssociation/DiseasePvalue.csv for plot '''

  ##########################  
  GO<-read.csv("DiseasePvalue.csv")
  head(GO)
  colnames(GO)<-c("CellType","DiseaseCat","pValue","Odds_ratio","JaccardIndex")
  library(ggplot2)
  #jpeg("ASD.Disease.jpeg",width = 800,height = 400)
  
ggplot(GO, aes(x=factor(CellType,
                               level=c("Astrocytes","OPC","OLIG","AST","MG")),
                      y=factor(DiseaseCat,
                               level=rev(c("ASD","AD","SCZ","EPI","BP","Anxiety","MDD")))))+
    theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    geom_point(aes(size=Odds_ratio, col=-log10(pValue)))+
    scale_size_continuous(range = c(1,5)) +
    scale_color_continuous(low = "grey", high = "red2")

ggsave("PvalueDisease.pdf",width=3.2,height=3)

range(-log10(GO$pValue))
range(GO$pValue)
  #dev.off()
  ##########################  ##########################  ##########################  ##########################    
##########################  ##########################  ##########################  ##########################  
##########################  ##########################  ##########################  ##########################  
##########################  ##########################  ##########################  ##########################   
