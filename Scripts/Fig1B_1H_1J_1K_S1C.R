
library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)
#----------------------------------------------------------------------------------------------------------
#### Fig.1B, Fig.1H, Fig.S1C Fig.1J, Fig.1K####
#----------------------------------------------------------------------------------------------------------
load("celltype_2.colors.rda")
load("Cell metaData.rda")
metaData$sample2 <- paste0(metaData$patient,"-",metaData$Field)
metaData$celltype_2 <- factor(metaData$celltype_2,levels = c("NaiveB","MemB-1","MemB-2","Plasmablast","Plasma-IgD",
                                                             "Plasma-IgM","Plasma-IgA1","Plasma-IgG1","Plasma-IgG3",
                                                             "Plasma-IgA2","Plasma-IgG2","Stressed plasma"))
# Fig.S1C
{
  df <- metaData %>% group_by(sample2,celltype_2,.drop = FALSE) %>% summarise(n=n())
  iOrd <-  paste0("P",1:16,"-"); iOrd2 <- NULL
  for (i in iOrd) {iOrd2 <- c(iOrd2,paste0(i,c("tumor","adjacent","intermediate","distant")))}
  df$sample2 <- factor(df$sample2,levels = iOrd2)
  ggplot(df, aes(x = sample2, y = n, fill = celltype_2)) + 
    geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
    scale_fill_manual("legend", values = celltype_2.colors) +
    theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Fig.1H
{
  df <- metaData %>% group_by(sample2,Field,patient,celltype_2,.drop = FALSE) %>% summarise(n=n())
  df <- rbind(df[df$Field=="tumor",],df[df$Field=="adjacent",],df[df$Field=="intermediate",],df[df$Field=="distant",])
  df$sample2 <- factor(df$sample2,levels = rev(unique(paste0(df$sample2))))
  ggplot(df, aes(x = sample2, y = n, fill = celltype_2)) + 
    geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
    scale_fill_manual("legend", values = celltype_2.colors) +
    theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

# Fig.1B
{
  df <- metaData %>% group_by(sample2,Field,.drop = FALSE) %>% summarise(n=n())
  iOrd <- unique(metaData[,c("sample2","EpiNegCells")])
  rownames(iOrd) <- iOrd$sample2
  df$EpiNegCells <- iOrd[paste0(df$sample2),"EpiNegCells"]
  df$frac <- df$n/df$EpiNegCells
  iOrd <-  paste0("P",1:16,"-"); iOrd2 <- NULL
  for (i in iOrd) {iOrd2 <- c(iOrd2,paste0(i,c("tumor","adjacent","intermediate","distant")))}
  df$sample2 <- factor(df$sample2,levels = iOrd2)
  my.colors <- c(tumor="#8c6bb1",adjacent="#4292c6",intermediate="#9ecae1",distant="#deebf7")
  ggplot(df, aes(x = sample2, y = frac, fill = Field)) + 
    geom_bar(position = "stack",stat = "identity") + #position="stack" gives numbers
    scale_fill_manual("legend", values = my.colors) +
    theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#Fig.1J
{
  df <- metaData[metaData$Field=="tumor",]
  df$celltype_2 <- factor(df$celltype_2,levels = c("NaiveB","MemB-1","MemB-2","Plasmablast","Plasma-IgD",
                                                               "Plasma-IgM","Plasma-IgA1","Plasma-IgG1","Plasma-IgG3",
                                                               "Plasma-IgA2","Plasma-IgG2","Stressed plasma"))
  df <- df %>% group_by(Mut,sample,celltype_2,.drop = FALSE) %>% summarise(n=n())
  iOrd <- unique(df$sample)# in the order of first occurrence in df.
  df$sample <- factor(df$sample,levels = iOrd)
  ggplot(df, aes(x = sample, y = n, fill = celltype_2)) + 
    geom_bar(position = "fill",stat = "identity") + #position="stack" gives numbers
    scale_fill_manual("legend", values = celltype_2.colors) +
    theme_classic()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  #EGFR mutated: P1-T2 P15T P6_T1 P7_T1 P8_T2
  #KRAS mutated: P10T P14T P2-T7
  #Others: P11T P12T P13T P16T P3-T6 P4-T2 P5-T2 P9_T1
}

#Fig.1K:
{
  df <- metaData[metaData$Field=="tumor",]
  df$celltype_1[df$celltype_2 %in% c("MemB-1","MemB-2")] <- "MemB"
  df$celltype_1[df$celltype_2 %in% c("NaiveB")] <- "NaiveB"
  df <- df %>% group_by(Mut,sample,celltype_1,.drop = FALSE) %>% summarise(n=n()) %>% mutate(freq = n/sum(n))
  iOrd <- unique(df$sample)# in the order of first occurrence in df.
  df$sample <- factor(df$sample,levels = iOrd)
  ggplot(df, aes(x = Mut, y = freq,fill = celltype_1)) + 
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(width = 0.2,size=2)+
    theme_classic()+ylab("% TME") + facet_wrap(~celltype_1)
}




