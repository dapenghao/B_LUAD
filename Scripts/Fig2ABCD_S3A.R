library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)

#----------------------------------------------------------------------------------------------------------
#### Fig.2A-D####
#----------------------------------------------------------------------------------------------------------
Bcell <- readRDS("Bcell.rds")
Idents(Bcell) <- Bcell$Field
load("celltype_2.colors.rda")
load("10Xclone_mutation.rda")
location.colors <- c(tumor="#8c6bb1",adjacent="#4292c6",intermediate="#9ecae1",distant="#deebf7")

# Fig.2D:
{
  library(CytoTRACE)
  Plasma <- subset(Bcell,subset=celltype_1=="Plasma")
  Plasma <- subset(Plasma,subset=celltype_2 != "Stressed plasma")
  GeneCounts <- as.matrix(Plasma@assays$RNA@counts)
  iOrd <- rowSums(GeneCounts>0)
  GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
  Plasma.CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE, ncores = 5, subsamplesize = 1000)
  Plasma$CytoTRACE_score <- Plasma.CytoTRACE$CytoTRACE[colnames(Plasma)]
  Plasma$celltype_2[Plasma$celltype_2 %in% c("Plasma-IgM","Plasma-IgD")] <- "Plasma-IgD"
  Plasma$IgH <- patient_obs[colnames(Plasma),"c_call"]
  df <- Plasma@meta.data[Plasma$IgH %ni% c("IGHE","IGHG4",NA),]
  df$IgH[df$IgH %in% c("IGHD","IGHM")] <- "IgM/D"
  df$IgH <- factor(df$IgH,levels=c("IgM/D","IGHG3","IGHG1","IGHA1","IGHG2","IGHA2"))
  p1 <- ggplot(df,aes(IgH, CytoTRACE_score)) + stat_compare_means()+
    geom_boxplot(aes(fill = IgH),outlier.shape = NA,alpha = 1)+
    theme_classic()+ylab("CytoTRACE score")
}

#1. Fig.2A-C and Fig.S3A
library(monocle3)
{
  PC <- subset(Bcell,subset=celltype_1 != "B")
  PC@active.assay <- "RNA"
  DEGs <- FindMarkers(PC,ident.1 = "Plasmablast",group.by = "celltype_2")
  DEGs$gene <- rownames(DEGs)
  DEGs <- DEGs[!grepl(pattern = "^IG[LHKJ]", DEGs$gene),]
  DEGs <- DEGs[!grepl(pattern = "^AC[0-9]", DEGs$gene),]
  DEGs <- DEGs[order(DEGs$p_val_adj),]
  cell_metadata <- PC@meta.data
  expression_data <- GetAssayData(PC,slot="counts")
  gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                                stringsAsFactors = F,row.names = rownames(expression_data))
  cds <- new_cell_data_set(expression_data=expression_data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = 30,norm_method="none",use_genes = DEGs$gene[1:500],method="PCA")
  cds <- reduce_dimension(cds,reduction_method="UMAP",preprocess_method="PCA",umap.min_dist = 0.35,cores=4,umap.n_neighbors = 15)
  cds <- cluster_cells(cds,partition_qval=0.05)
  cds <- learn_graph(cds,use_partition=F)
  cds <- order_cells(cds)
  
  # plot the expression for "STAT3","IKZF3","CD37","IGHA1"
  iOrd <- paste0(cds@clusters$UMAP$clusters)
  PC$pseudoTime <- cds@principal_graph_aux$UMAP$pseudotime
  PC$BCR <- colnames(PC) %in% rownames(patient_obs)
  cds@colData$BCR <- PC$BCR
  
  # CytoTrace Score:
  {
    library(CytoTRACE)
    GeneCounts <- as.matrix(PC@assays$RNA@counts)
    iOrd <- rowSums(GeneCounts>0)
    GeneCounts <- GeneCounts[iOrd>10,]#only keep genes expressing in more than 10 cell
    PC.CytoTRACE <- CytoTRACE(GeneCounts, enableFast = TRUE, ncores = 5, subsamplesize = 1000)# run in server
  }
  PC.metaData <- PC@meta.data
  PC.metaData$CytoTrace <- PC.CytoTRACE$CytoTRACE[rownames(PC.metaData)]
  
  # plot gene expression across trajectory:
  {
    genes <- c("STAT3","IKZF3","SDC1","IGHA1","MS4A1","TYMS","MKI67","TUBB","CD37","IGHG1","IGHG2","IGHA2")
    for (i in genes) {
      print(plot_cells(cds, label_groups_by_cluster=F, genes=i,
                       reduction_method="UMAP",show_trajectory_graph=F,
                       label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
                       alpha=1,cell_size = 1,group_label_size = 4))#+scale_colour_continuous(low="#F7FCFD",high="#8C6BB1")
    }
  }
  
  # statistics along PseudoTime:
  {
    # Smooth data for fraction of cells with BCR posotivity across trajectory:
    df <- PC.metaData
    df <- df[order(df$pseudoTime),]
    ggplot(df, aes(x=pseudoTime,y=as.numeric(BCR)))+geom_smooth(fullrange=T,method="loess",span=0.5)+
      labs(title="TRUE") +xlim(0,51)+ theme_classic()
    
    # Smooth data for fraction of cells with Smoking History across trajectory:
    df <- PC.metaData
    df <- df[order(df$pseudoTime),]
    ggplot(df, aes(x=pseudoTime,y=as.numeric(Smoking2=="smoker")))+geom_smooth(fullrange=T,method="loess",span=0.5)+
      labs(title="TRUE") +xlim(0,51)+ theme_classic()
    
    # Smooth data for CytoTrace Score:
    df <- PC.metaData
    df <- df[order(df$pseudoTime),]
    ggplot(df, aes(x=pseudoTime,y=CytoTrace))+geom_smooth(fullrange=T,method="loess",span=1)+
      labs(title="TRUE") +xlim(0,50)+ theme_classic()
    cds@colData$CytoTrace <- PC.metaData$CytoTrace
    plot_cells(cds, label_groups_by_cluster=F, color_cells_by = "CytoTrace",
               reduction_method="UMAP",show_trajectory_graph=T,alpha=1,
               label_roots = F,label_leaves = F,label_cell_groups =F,label_branch_points = F,
               cell_size = 0.5,group_label_size = 4)
  }
}
