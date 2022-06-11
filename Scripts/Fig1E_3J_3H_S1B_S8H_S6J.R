library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)

Bcell <- readRDS("Bcell.rds")
load("celltype_2.colors.rda")
#----------------------------------------------------------------------------------------------------------
#### Fig.1E, Fig.S1B ####
#----------------------------------------------------------------------------------------------------------
load("10Xclone_mutation.rda")
# Fig.1E
{
  p1 <- DimPlot(Bcell, reduction = "umap",label =F,group.by = "celltype_2",cols=celltype_2.colors,pt.size=1)
  p1 <- AugmentPlot(p1)
  
  my.colors <- c(smoker="#cb181d",`never smoker`="#fee0d2")
  p2 <- DimPlot(Bcell, reduction = "umap",label =F,group.by = "Smoking2",cols=my.colors,pt.size=1,shuffle=T)
  p2 <- AugmentPlot(p2)
  
  my.colors <- c(tumor="#8c6bb1",adjacent="#4292c6",intermediate="#9ecae1",distant="#deebf7")
  p3 <- DimPlot(Bcell, reduction = "umap",label =F,group.by = "Field",cols=my.colors,pt.size=1,shuffle=T)
  p3 <- AugmentPlot(p3)
  
  
  load("./10Xclone_mutation.rda")
  iOrd <- intersect(patient_obs$cell_id,colnames(Bcell))
  Bcell.BCR <- subset(Bcell,cells=iOrd)
  Bcell.BCR$Isotype <- patient_obs[colnames(Bcell.BCR),"c_call"]
  Bcell.BCR$MF_group <- patient_obs[colnames(Bcell.BCR),"MF_group"]
  Bcell.BCR$Field <- ifelse(Bcell.BCR$Field=="tumor","tumor","normal")
  Bcell.BCR$Mut <- paste0(Bcell.BCR$Field,"-",Bcell.BCR$Mut)
  Bcell.noBCR <- subset(Bcell,cells=setdiff(colnames(Bcell),patient_obs$cell_id))
  Bcell.noBCR$Isotype="noBCR"
  
  isotype.colors <- c(IGHA1="#fec44f",IGHA2="#fe9929",IGHG1="#c6dbef",IGHG2="#6baed6",IGHG3="#2171b5",IGHG4="#08306b",
                      IGHD="#807dba",IGHE="#D51F26",noBCR="#d9d9d9",IGHM="#fa9fb5")
  p4 <- DimPlot(Bcell.BCR, reduction = "umap",label =F,group.by = "Isotype",cols=isotype.colors,pt.size=1.5)
  p4 <- AugmentPlot(p4)
}

# Fig.S1B
{
  p1 <- FeaturePlot(Bcell,features = c("FCER2","MZB1","JCHAIN","AIM2"))
  p1 <- AugmentPlot(p1)
}
#----------------------------------------------------------------------------------------------------------
#### Fig.3J ####
#----------------------------------------------------------------------------------------------------------
load("10Xclone_mutation.rda")
patient_obs$MF_group[patient_obs$MF_group %in% c("Medium","High")] <- "High"
ggplot(patient_obs, aes(x = MF_group, y = log10(clone_size),fill = MF_group)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values=my.colors)+
  theme_classic()+ylab("Clone size")+ylim(c(0,3.5))
#----------------------------------------------------------------------------------------------------------
#### Fig.3H, Fig.S8H ####
#----------------------------------------------------------------------------------------------------------
load("10Xclone_mutation.rda")
#Fig.3H:
{
  iOrd <- intersect(patient_obs$cell_id,colnames(Bcell))
  Bcell.mut <- subset(Bcell,cells=iOrd)
  Bcell.mut$MF <- patient_obs[colnames(Bcell.mut),"mu_freq_seq_r"]
  Bcell.mut$MF_group <- Bcell.mut$MF
  Bcell.mut$MF_group <- case_when(
    Bcell.mut$MF == 0 ~ "Germline",
    Bcell.mut$MF <=0.03 ~ "Low",
    #Bcell.mut$MF <=0.06 ~ "Medium",
    Bcell.mut$MF >0.03  ~ "High")
  
  my.colors <- c(Germline="#74add1",Low="#ffffbf",Medium="#fdae61",High="#d73027")
  p1 <- DimPlot(Bcell.mut, reduction = "umap",label =F,group.by = "MF_group",cols=my.colors,pt.size=1)
  p1 <- AugmentPlot(p1)
}

#Fig.S8H:
{
  iOrd <- intersect(patient_obs$cell_id,colnames(Bcell))
  Bcell.mut <- subset(Bcell,cells=iOrd)
  Bcell.mut$MF <- patient_obs[colnames(Bcell.mut),"mu_freq_seq_r"]
  Bcell.mut$MF_group <- Bcell.mut$MF
  Bcell.mut$MF_group <- case_when(
    Bcell.mut$MF == 0 ~ "Germline",
    Bcell.mut$MF <=0.03 ~ "Low",
    Bcell.mut$MF >0.03  ~ "High")
  
  my.colors <- c(Germline="#74add1",Low="#ffffbf",High="#d73027")
  p1 <- DimPlot(Bcell.mut, reduction = "umap",label =F,group.by = "MF_group",cols=my.colors,pt.size=0.5,split.by = "Mut")
  p1 <- AugmentPlot(p1)
  pdf("png_UMAPs of cells by mutation groups splitted.pdf",width = 18,height=6)
  print(p1)
  dev.off()
}

#----------------------------------------------------------------------------------------------------------
#### Fig.S8J ####
#----------------------------------------------------------------------------------------------------------
library(scales)
load("10Xclone_mutation.rda")
my.colors <- c(Germline="#74add1",Low="#ffffbf",High="#d73027")
patient_obs <- patient_obs[patient_obs$Field=="tumor",]
patient_obs$MF_group[patient_obs$MF_group %in% c("Medium","High")] <- "High"
x <- table(patient_obs$MutGenes)/nrow(patient_obs)
expected <- as.numeric(x);names(expected) <- names(x)
table(patient_obs$MF_group)
df <- NULL
MF_group <- unique(paste0(patient_obs$MF_group))
for (i in MF_group) {
  x <- patient_obs[patient_obs$MF_group==i,]
  x <- table(x$MutGenes)/nrow(x)/expected
  x <- data.frame(MF_group=i,MutGenes=names(x),Roe=as.numeric(x))
  df <- rbind(df,x)
}

labs <- as.character(round(df$Roe,2))
df2 <- df;df2$Roe[df2$Roe>=2.5]=2.5#change the max value to 2
df2$MF_group <- factor(df2$MF_group,levels=names(my.colors))
p1 <- ggplot(data =  df2, aes(x = MutGenes, y = MF_group)) + 
  geom_tile(aes(fill = Roe)) +
  scale_fill_gradientn(colours=c("#41ae76","#99d8c9","#f0f0f0","#feb24c","#fc4e2a"),
                       values=rescale(c(0.5,0.75,1,1.75,2.5)),
                       guide="colorbar")
p1 +  geom_text(aes(label=labs), size=4) + theme_classic()



