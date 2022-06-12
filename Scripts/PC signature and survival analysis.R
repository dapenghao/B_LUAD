
library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)
Bcell <- readRDS("./Bcell.rds")
Bcell$celltype_1[Bcell$celltype_2 %in% c("MemB-1","MemB-2")] <- "MemB"
Bcell$celltype_1[Bcell$celltype_2 %in% c("NaiveB")] <- "NaiveB"
Bcell@active.assay <- "RNA"
#----------------------------------------------------------------------------------------------------------
#### GSE93157 (Prat A. et.al. dataset) ####
#----------------------------------------------------------------------------------------------------------

#1. processing the data:
{
  library(GEOquery)
  gse <-getGEO("GSE93157")
  meta <- pData(gse[[1]])
  meta <- meta[grepl("LUNG",meta$source_name_ch1),]
  meta$`best.resp:ch1`[meta$`best.resp:ch1` %in% c("CR","PR")] <- "CRPR"
  
  exp.mat <- exprs(gse[[1]])
  exp.mat <- exp.mat[,rownames(meta)]
  iOrd <- rowSums(is.na(exp.mat))
  exp.mat <- exp.mat[-1*which(iOrd==35),]
}

#2. MCPcounter estimate:
{
  library(MCPcounter)
  PC.DEGs <- c("JCHAIN","TNFRSF17","MZB1","DERL3","SPATS2","JSRP1","TXNDC5","MANEA")
  B.DEGs <- data.frame(DEGs = PC.DEGs,cellType = rep("PC",length(PC.DEGs)))
  B.DEGs <- data.frame(DEGs = c(MemB.DEGs,NaiveB.DEGs,PC.DEGs),
                       cellType = c(rep("MemB",length(MemB.DEGs)),
                                    rep("NaiveB",length(NaiveB.DEGs)),
                                    rep("PC",length(PC.DEGs))))
  MCPEstimates=MCPcounter.estimate(exp.mat,featuresType="affy133P2_probesets",probesets=B.DEGs)
  MCPEstimates <- t(MCPEstimates)
  meta$PC.MCP <- MCPEstimates[,'PC']
  meta$PC.MCP <- scale(meta$PC.MCP)[,1]
}

#3. survival analysis
{
  library(survival);library(survminer)
  meta$Group <- ifelse(meta$`best.resp:ch1` %in% c("CRPR","SD"),"DCB","NCB")
  meta$Group <- NA
  meta$Group[(meta$`best.resp:ch1` %in% c("CRPR","SD")) & meta$PFS > 24*7] <- "DCB"
  meta$Group[(meta$`best.resp:ch1` %in% c("PD","SD")) & meta$PFS <= 24*7 & meta$Recurrence==1] <- "NDB"
  #meta$Group[(meta$`best.resp:ch1` %in% c("PD","SD")) & meta$PFS <= 24*7] <- "NDB"
  
  meta$PFS <- gsub('pfs: ','',meta$characteristics_ch1.14)
  meta$PFS <- as.numeric(meta$PFS)
  meta$Recurrence <- gsub("pfse: ","",meta$characteristics_ch1.13)
  meta$Recurrence <- as.numeric(meta$Recurrence)
  meta$PFS <- meta$PFS*30.5
  
  fit <- survfit(Surv(PFS, Recurrence) ~ meta$PC.MCP>0,data = meta)
  ggsurvplot(fit, data = meta, risk.table = TRUE,pval = TRUE,risk.table.y.text = FALSE,
             xscale=30.5,break.x.by=30.5*6,xlim=c(0,30.5*25),
             palette = c("#2166ac","#d6604d"))
}

#----------------------------------------------------------------------------------------------------------
#### PC signature across stages and survival analysis of TCGA dataset ####
#----------------------------------------------------------------------------------------------------------
load("RNAseq data of LUADfrom xena.rda") 
load("Clin of LUAD.rda") # TCGA data were obtained from Xena Browser, see Supplementary Methods for details.
Clin.LUAD <- Clin.LUAD[colnames(RNAseq.xena),]# only keep samples with RNAseq
RNAseq.xena <- RNAseq.xena[,rownames(Clin.LUAD)]

#1. MCPcounter estimates of PC signature:
PC.DEGs <- c("JCHAIN","TNFRSF17","MZB1","DERL3","SPATS2","JSRP1","TXNDC5","MANEA")
B.DEGs <- data.frame(DEGs = PC.DEGs,cellType = rep("PC",length(PC.DEGs)))
B.DEGs$DEGs[B.DEGs$DEGs=="JCHAIN"] <- "IGJ"
MCPEstimates=MCPcounter.estimate(RNAseq.xena,featuresType="affy133P2_probesets",probesets=B.DEGs)
MCPEstimates <- t(MCPEstimates)
Clin.LUAD$PC.MCP <- MCPEstimates[,'PC']
Clin.LUAD$PC.MCP <- scale(Clin.LUAD$PC.MCP)[,1]

#2. plot across stages:
Clin.LUAD$Stage <- case_when(
  Clin.LUAD$tumor_stage.diagnoses %in% c("stage i","stage ia","stage ib") ~ "I",
  Clin.LUAD$tumor_stage.diagnoses %in% c("stage ii","stage iia","stage iib") ~ "II",
  Clin.LUAD$tumor_stage.diagnoses %in% c("stage iiia","stage iiib") ~ "III",
  Clin.LUAD$tumor_stage.diagnoses %in% c("stage iv") ~ "IV",
  Clin.LUAD$tumor_stage.diagnoses %in% c("not reported") ~ "NA",
)

ggplot(Clin.LUAD, aes(x = Stage, y = PC.MCP,fill = Stage)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2,size=2)+
  theme_classic()+ylab("PC signature")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(legend.position = "none")

#3. smoothHR analysis:
library(phenoTest)
smoothCoxph(Clin.LUAD$OS.time,Clin.LUAD$OS,Clin.LUAD$PC.MCP,xlab="PC signature",main="PC signature")
gene.cox <- coxph(Surv(OS.time, OS==1) ~ PC.MCP, Clin.LUAD)
summary(gene.cox)

#4. survival analysis:
df <- Clin.LUAD
df$Group <- ifelse(df$PC.MCP>0.5,"high","middle")
df$Group[df$PC.MCP < -0.5] <- "low"
fit <- survfit(Surv(OS.time, OS==1) ~ Group,data = df)
ggsurvplot(fit, data = df, risk.table = TRUE,pval = TRUE,risk.table.y.text = FALSE,
           xscale=365,break.x.by=2*365,xlim=c(0,365*6),
           palette = c("#ef6548","#0570b0","#bdbdbd"))
survdiff(Surv(OS.time, OS==1) ~ Group,data = df[df$Group!="middle",])

