setwd("/Users/dhao/Desktop/projects/Bcell_lung")
library(ggplot2);library(data.table);library(Seurat);
library(dplyr);library(tidyr);library(ggpubr);library(RColorBrewer)

Bcell <- readRDS("Bcell.rds")
load("celltype_2.colors.rda")
Bcell$celltype_2 <- factor(Bcell$celltype_2,levels = c("NaiveB","MemB-1","MemB-2","Plasmablast","Plasma-IgD",
                                                       "Plasma-IgM","Plasma-IgA1","Plasma-IgG1","Plasma-IgG3",
                                                       "Plasma-IgA2","Plasma-IgG2","Stressed plasma"))
#----------------------------------------------------------------------------------------------------------
#### Fig.1F ####
#----------------------------------------------------------------------------------------------------------
x=c("MS4A1","CD19","BANK1","SELL","IL4R","FCER2","TCL1A","BACH2",#Pan-B and Naive B cells
    "CD27","CD38","JCHAIN","MZB1","SDC1","TNFRSF17",#Plasma and memory B (CD27)
    "IGHM", "IGHD","IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHA1", "IGHA2",#plasma isotype
    "MKI67","TUBB","STMN1","TYMS",# plasmablast
    "BAG3","HSPA6","HSPB1"# Stressed Plasma
)
DotPlot(Bcell, features = x,assay = "RNA",scale = T,group.by = "celltype_2") + 
  scale_colour_gradientn(colors=brewer.pal(10, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=x,labels=x)+ ylab("thres>0")
#----------------------------------------------------------------------------------------------------------
#### Fig.S1A ####
#----------------------------------------------------------------------------------------------------------
MHC2 <- c("CD74","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1",
          "HLA-DPB1","HLA-DQA1","HLA-DQA2","HLA-DQB1","HLA-DQB2",
          "HLA-DRA","HLA-DRB1","HLA-DRB3","HLA-DRB4","HLA-DRB5")
MHC2 <- MHC2[MHC2 %in% rownames(Bcell)]
DotPlot(object = Bcell, features = MHC2,scale = T,group.by = "celltype_2") + 
  scale_colour_gradientn(colors=brewer.pal(10, "BuPu")) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=MHC2,labels=MHC2)#MHC2
#----------------------------------------------------------------------------------------------------------
#### Fig.11B ####
#----------------------------------------------------------------------------------------------------------
load("./Revision/Github/Input_data/MarkerExp.rda")
load("./Revision/Github/Input_data/cellType_All.rda")
df <- unique(Bcell@meta.data[,c("orig.ident","Field")]);rownames(df) <- df$orig.ident
cellType_All$Field <- df[cellType_All$orig.ident,"Field"]
obj <- CreateSeuratObject(t(MarkerExp),assay = "RNA",meta.data = cellType_All)
obj@assays$RNA@data <- as.matrix(t(MarkerExp))
obj <- subset(obj,subset=cellType %in% c("CD4T","CD8T","Fibroblast","Myeloid"))
obj$cellType <- paste0(obj$Field,"-",obj$cellType)
p1 <- DotPlot(object = obj, features = colnames(MarkerExp),assay = "RNA",scale = T,group.by = "cellType") + 
  scale_colour_gradientn(colors=brewer.pal(9, "BuPu"))  + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=colnames(MarkerExp),labels=colnames(MarkerExp))

#----------------------------------------------------------------------------------------------------------
#### Fig.11D ####
#----------------------------------------------------------------------------------------------------------
Bcell$celltype_1[Bcell$celltype_2 %in% c("MemB-1","MemB-2")] <- "MemB"
Bcell$celltype_1[Bcell$celltype_2 %in% c("NaiveB")] <- "NaiveB"
Bcell$celltype_1[grep("IgA",Bcell$celltype_2)] <- "IgA-Plasma"
Bcell$celltype_1[grep("IgG",Bcell$celltype_2)] <- "IgG-Plasma"
Bcell2 <- subset(Bcell,subset=celltype_1 %in% c("NaiveB","MemB","IgA-Plasma","IgG-Plasma"))
Bcell2$celltype_1 <- factor(Bcell2$celltype_1,levels = c("NaiveB","MemB","IgA-Plasma","IgG-Plasma"))
sel.genes <- c("CCR7","P2RY8","GPR183","TGFB1","CD24","NR4A1",
               "NR4A2","LTB","BTLA","CCR10","VSIR","TNFRSF18","LGALS3")
DotPlot(object = Bcell2, features = sel.genes,assay = "RNA",scale = T,group.by = "celltype_1") + 
  scale_colour_gradientn(colors=brewer.pal(9, "BuPu"))  + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(breaks=sel.genes,labels=sel.genes)



