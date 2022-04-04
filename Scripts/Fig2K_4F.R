library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)

#----------------------------------------------------------------------------------------------------------
#### Fig.2K####
#----------------------------------------------------------------------------------------------------------
Bcell <- readRDS("Bcell.rds")
MemB <- subset(Bcell,subset=celltype_2 %in% c("MemB-1","MemB-2"))
iOrd <- c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3")
iOrd <- GetAssayData(MemB,slot = "data")[iOrd,]
MemB$IgAG <- colMeans(iOrd)
MemB$development <- "NA"
MemB$development[MemB$IgAG>=0 & MemB$IgAG<1] <- "Early"
MemB$development[MemB$IgAG>=1 & MemB$IgAG<2] <- "Intermediate"
MemB$development[MemB$IgAG>=2] <- "Late"

DEGs_LvE <- FindMarkers(MemB,ident.1="Late",ident.2 = "Early",group.by = "development",logfc.threshold =0)
MemB$Field2 <- ifelse(MemB$Field=="tumor","tumor","normal")
DEGs_TvN <- FindMarkers(MemB,ident.1="tumor",ident.2 = "normal",group.by = "Field2",logfc.threshold =0)

MemB2<- subset(MemB,cells=intersect(colnames(MemB),rownames(patient_obs)))
MemB2$group <- patient_obs[colnames(MemB2),"c_call"]
MemB2$group[MemB2$group %in% c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4")] <- "Switched"
MemB2$group[MemB2$group %in% c("IGHM","IGHD")] <- "Unswitched"
DEGs_SvU <- FindMarkers(MemB2,ident.1="Switched",ident.2 = "Unswitched",group.by = "group",logfc.threshold =0)

colnames(DEGs_SvU) <- paste0(colnames(DEGs_SvU),"-","SvU")
colnames(DEGs_TvN) <- paste0(colnames(DEGs_TvN),"-","TvN")
colnames(DEGs_LvE) <- paste0(colnames(DEGs_LvE),"-","LvE")

iOrd <- intersect(rownames(DEGs_SvU),rownames(DEGs_LvE))
df <- cbind(DEGs_SvU[iOrd,],DEGs_LvE[iOrd,])
df$label <- rownames(df)
df$label[abs(df$`avg_log2FC-SvU`)<0.5 | abs(df$`avg_log2FC-LvE`)<0.5] <- NA
df$log2FC_TvN <- DEGs_TvN[iOrd,"avg_log2FC-TvN"]
df$log2FC_TvN[df$log2FC_TvN> 1] <- 1
df$log2FC_TvN[df$log2FC_TvN< -1] <- -1

library(ggrepel)
ggplot(data=df, aes(x = `avg_log2FC-SvU`, y = `avg_log2FC-LvE`, col=log2FC_TvN,label=label)) +
  geom_point(stroke = 0, alpha = 0.8,shape = 16,size=3) + 
  theme_classic() +
  geom_text_repel(colour="black") +
  scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")

#----------------------------------------------------------------------------------------------------------
#### Fig.4F####
#----------------------------------------------------------------------------------------------------------
Bcell <- readRDS("Bcell.rds")
Bcell$Field2 <- ifelse(Bcell$Field=="tumor","tumor","normal")
MemB <- subset(Bcell,subset=celltype_2 %in% c("MemB-1","MemB-2"))
MemB <- subset(Bcell,subset=celltype_2 %in% c("Plasma"))
NaiveB <- subset(Bcell,subset=celltype_1=="NaiveB")

DEGs_TvN <- FindMarkers(MemB,ident.1="tumor",ident.2 = "normal",group.by = "Field2",logfc.threshold =0)
DEGs_TvN_MemB <- DEGs_TvN
DEGs_TvN_PC <- FindMarkers(Plasma,ident.1="tumor",ident.2 = "normal",group.by = "Field2",logfc.threshold =0)
DEGs_TvN_Naive <- FindMarkers(NaiveB,ident.1="tumor",ident.2 = "normal",group.by = "Field2",logfc.threshold =0)

colnames(DEGs_TvN_Naive) <- paste0(colnames(DEGs_TvN_Naive),"_","Naive")
colnames(DEGs_TvN_MemB) <- paste0(colnames(DEGs_TvN_MemB),"_","MemB")
colnames(DEGs_TvN_PC) <- paste0(colnames(DEGs_TvN_PC),"_","PC")

iOrd <- intersect(rownames(DEGs_TvN_MemB),rownames(DEGs_TvN_PC))
iOrd <- intersect(iOrd,rownames(DEGs_TvN_Naive))

df <- cbind(DEGs_TvN_MemB[iOrd,],DEGs_TvN_Naive[iOrd,])
df$label <- rownames(df)
df$label[abs(df$avg_log2FC_MemB)<0.5 | abs(df$avg_log2FC_Naive)<0.5] <- NA
df$log2FC_PC <- DEGs_TvN_PC[iOrd,"avg_log2FC_PC"]
df$log2FC_PC[df$log2FC_PC> 1] <- 1
df$log2FC_PC[df$log2FC_PC< -1] <- -1

library(ggrepel)
ggplot(data=df, aes(x = `avg_log2FC_MemB`, y = `avg_log2FC_Naive`, col=log2FC_PC,label=label)) +
  geom_point(stroke = 0, alpha = 0.8,shape = 16,size=3) + 
  theme_classic() +
  geom_text_repel(colour="black") +
  scale_color_gradient2(midpoint=0, low="#313695", mid="#f0f0f0",high="#a50026",space ="Lab") +
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  geom_vline(xintercept=0, linetype="dashed", color = "red")




