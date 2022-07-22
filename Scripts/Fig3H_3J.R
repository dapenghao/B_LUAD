library(ggplot2);library(data.table);library(Seurat);library(dplyr);library(tidyr);library(ggpubr)

#----------------------------------------------------------------------------------------------------------
#### Fig.3H####
#----------------------------------------------------------------------------------------------------------
Bcell <- readRDS("Bcell.rds")
Idents(Bcell) <- Bcell$Field
MemB <- subset(Bcell,subset=celltype_2 %in% c("MemB-1","MemB-2"))
#monocle2:
{
  library(monocle)
  iOrd <- c("IGHA1","IGHA2","IGHG1","IGHG2","IGHG3")
  iOrd <- GetAssayData(MemB,slot = "data")[iOrd,]
  MemB$IgAG <- colMeans(iOrd)
  MemB$development <- "NA"
  MemB$development[MemB$IgAG>=0 & MemB$IgAG<1] <- "Early"
  MemB$development[MemB$IgAG>=1 & MemB$IgAG<2] <- "Intermediate"
  MemB$development[MemB$IgAG>=2] <- "Late"
  
  cell_metadata <- MemB@meta.data
  expression_data <- GetAssayData(MemB,slot="counts")
  gene_annotation <- data.frame(gene_short_name=rownames(expression_data),
                                stringsAsFactors = F,row.names = rownames(expression_data))
  pd <- new("AnnotatedDataFrame", data = cell_metadata)
  fd <- new("AnnotatedDataFrame", data = gene_annotation)
  MemB.cds <- newCellDataSet(cellData=expression_data,phenoData = pd,
                             featureData = fd,expressionFamily=negbinomial.size())
  expression_data <- expression_data[rowSums(expression_data>0)>10,]#only keep genes expressing in more than 10 cell
  MemB.cds$Field2 <- MemB.cds$Field
  MemB.cds$Field2[MemB.cds$Field2 %in% c("adjacent","distant","intermediate")] <- "normal"
  
  ####### Run monocle ##############################
  MemB.cds <- estimateSizeFactors(MemB.cds)
  MemB.cds <- estimateDispersions(MemB.cds)
  MemB.cds <- detectGenes(MemB.cds, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(MemB.cds),num_cells_expressed >= 300))
  diff_test_res <- differentialGeneTest(MemB.cds[expressed_genes,],fullModelFormulaStr = "~development")
  ordering_genes <- unique(row.names(diff_test_res)[order(diff_test_res$qval)][1:1000])
  MemB.cds <- setOrderingFilter(MemB.cds, ordering_genes)
  MemB.cds <- reduceDimension(MemB.cds, max_components = 2,reduction_method = 'DDRTree',pseudo_expr=1)
  MemB.cds <- orderCells(MemB.cds)
  
  plot_cell_trajectory(MemB.cds, color_by = "Field2",cell_size=0.1)
  plot_cell_trajectory(MemB.cds, color_by = "CytoTRACE_score",cell_size=0.5)
  plot_cell_trajectory(MemB.cds, color_by = "IgM",cell_size=0.5)
  plot_cell_trajectory(MemB.cds, color_by = "IgAG",cell_size=0.5)
  plot_cell_trajectory(MemB.cds, color_by = "celltype_2",cell_size=0.5)
  plot_cell_trajectory(MemB.cds, color_by = "development",cell_size=0.5)
  ####### end ##############################
}

#----------------------------------------------------------------------------------------------------------
#### Fig.3J####
#----------------------------------------------------------------------------------------------------------
library(smoother);library(pheatmap);library(Hmisc)
df=data.frame(t(MemB.cds@reducedDimS))
df$Pseudotime <- MemB.cds$Pseudotime
df <- df[order(df$Pseudotime),]
df$group <- as.numeric(cut2(df$Pseudotime, g=100))

cor.mat <- GetAssayData(MemB,slot = "data")[,rownames(df)]
cor.mat <- as.matrix(cor.mat)
cor.mat <- groupMeans(cor.mat,groups = df$group)

# genes variable across pseudotime: -- use this one
{
  iOrd <- apply(cor.mat, 1, sd)
  cor.mat <- cor.mat[order(iOrd,decreasing = T),]
  result.x <- cor.mat[1:3000,]
  gene.cor <- data.frame(genes=rownames(result.x),cors=cor(1:100,t(result.x),use="na.or.complete")[1,])
  gene.cor <- gene.cor[order(gene.cor$cors),]
  gene.cor$BCRgenes <- NA
  gene.cor$BCRgenes[grep("^IG",gene.cor$genes)] <- "BCRgenes"
}

result.x.smth <- result.x[gene.cor$genes,]
for (i in 1:nrow(result.x)) {
  result.x.smth[i,] <- smth(result.x.smth[i,],window = 0.1,method = "gaussian")
}
result.x.smth[,1:5] <- result.x.smth[,6]
result.x.smth[,96:100] <- result.x.smth[,95]
bk <- c(seq(-1,-0.1,by=0.02),seq(0,1,by=0.02))
colour_bk <- c(colorRampPalette(c("#2166ac","#d1e5f0"))(38),
               colorRampPalette(c("#d1e5f0","#f7f7f7"))(10),
               colorRampPalette(c("#f7f7f7","#fddbc7"))(10),
               colorRampPalette(c("#fddbc7","#b2182b"))(39))
pheatmap(result.x.smth, cluster_rows=F, cluster_cols=F,
         scale="row",show_colnames=F,show_rownames=F,
         breaks = bk,color = colour_bk,angle_col=45,
         clustering_method="ward.D",border_color=NA,
         annotation_row=gene.cor[,"BCRgenes",drop=F])











