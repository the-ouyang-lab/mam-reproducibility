rm(list = ls())
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(cluster)
library(msigdbr)
library(ComplexHeatmap)
if(!dir.exists("res_1inteData/")){dir.create("res_1inteData/")}
hmapCol <- colorRampPalette(c("blue", "white", "red"))(21)
hmapBrk <- seq(-1, 1, length.out = 21)
genes.to.label <- c(
  'SELENOP','C1QA','RNASE1' ,'C1QC','FCN1', 'SDC2', 'PLA2G7','S100A9','S100A8',
  'VCAN', 'APOE','GPNMB','CD9','STAB1','SPP1', 'FABP5', 'TREM2', 'F13A1', 
  'CCL20', 'FOLR2', 'INHBA',  'CHI3L1',  'MATK','PTGS2','IL1B', 'FABP4', 'LYVE1', 
  'APOC1', 'CCL4','S100A12', 'CTSK', 'LST1', 'CDKN1C', 'FCGR3A',
  'CCL3', 'DAB2','CCL4L2', 'PHLDA3', 'EGR1', 'ATF3', 'VEGFA', 'SEMA6B', 
  'IFI44L', 'IFIT3', 'IFIT1', 'MT1G', 'MT1X', 'MT1H')



inpFiles = list.files("res_0prepData2rds/", pattern = ".rds")
for(i in inpFiles){
  tmp = readRDS(paste0("res_0prepData2rds/", i))
  tmpG = unique(tstrsplit(colnames(tmp), "_")[[1]])
  tmpG = c(length(tmpG), length(grep("norm", tmpG)), length(grep("norm", tmpG, invert = T)))
  print(paste0(i, ": ", tmpG[1], " patients | ", tmpG[2], " | ", tmpG[3]))
}

oupSil = data.table()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  ### 1. IO and preproc
  inpFiles = list.files("res_0prepData2rds/", pattern = ".rds")
  inpFiles = inpFiles[grep(i, inpFiles)]
  inpGEX = list()
  commonGenes = c()
  for(j in inpFiles){
    tmp = readRDS(paste0("res_0prepData2rds/", j))
    tmp = tmp[rowSums(tmp > 0) > 0, ] # get expressed genes
    commonGenes = c(commonGenes, rownames(tmp))
    inpGEX[[strsplit(j, "\\.")[[1]][1]]] = tmp
  }
  
  # Get genes that are expressed in 2 or more datasets
  commonGenes = names(table(commonGenes))[table(commonGenes) >= 2]
  for(j in names(inpGEX)){
    tmpNewG = setdiff(commonGenes, rownames(inpGEX[[j]]))
    if(length(tmpNewG) > 0){
      tmpNewM = matrix(0, ncol = ncol(inpGEX[[j]]), nrow = length(tmpNewG))
      rownames(tmpNewM) = tmpNewG; colnames(tmpNewM) = colnames(inpGEX[[j]])
      inpGEX[[j]] = rbind(inpGEX[[j]], tmpNewM)
    }
    inpGEX[[j]] = inpGEX[[j]][commonGenes, ]
  }
  
  # For heart, split heart2-norm
  if(i == "heart"){
    inpGEX[["oupGEXr3_heart2norm"]] = inpGEX[["oupGEXr3_heart2"]][
      , grep("heart2-norm", colnames(inpGEX[["oupGEXr3_heart2"]]))]
    inpGEX[["oupGEXr3_heart2"]] = inpGEX[["oupGEXr3_heart2"]][
      , grep("heart2-norm", colnames(inpGEX[["oupGEXr3_heart2"]]), invert = TRUE)]
  }
  
  
  ### 2. SCT and Seurat integration
  # Start SCTransform
  oupSeuList = list()
  for(j in names(inpGEX)){
    tmp = inpGEX[[j]]
    if(j == "oupGEXr3_kidney2"){names(colnames(tmp)) = NULL}
    tmp = CreateSeuratObject(counts = tmp, project = j)
    if(j == "oupGEXr3_heart2norm"){
      tmp = SCTransform(tmp, verbose = F)
    } else {
      tmp = SCTransform(tmp, vars.to.regress = "orig.ident", verbose = F)
    }
    tmp$study = j
    oupSeuList[[j]] = tmp
  }
  
  # Start integration
  features <- SelectIntegrationFeatures(object.list = oupSeuList, nfeatures = 3000)
  seu <- PrepSCTIntegration(object.list = oupSeuList, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = seu, normalization.method = "SCT",
                                    anchor.features = features)
  seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
  seu <- FindClusters(seu, resolution = (4:20)/20, verbose = FALSE)

  # Plot 
  colIdent = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$orig.ident))
  ggOut = list()
  ggOut[["sample"]] = DimPlot(seu, group.by = "orig.ident", cols = colIdent, 
                              shuffle = T, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tmpG = c("SPP1", "MMP7", "MMP9", "TREM2", "FCN1")
  for(j in tmpG){
    if(j %in% rownames(seu)){DefaultAssay(seu) = "integrated"}
    if(!j %in% rownames(seu)){DefaultAssay(seu) = "RNA"}
    ggOut[[j]] = FeaturePlot(seu, features = j, order = T, pt.size = 0.25) + 
      scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  ggsave(wrap_plots(ggOut), height = 12, width = 16,
         filename = paste0("res_1inteData/", i, "_umap.png"))
  
  
  ### 3. Merge and annotate clusters
  # Calculate silhouette score and determine best res
  if(i == "lung"){
    set.seed(42); idx = sample.int(ncol(seu), size = 25000)
  } else {
    idx = 1:ncol(seu)
  }
  tmpD = dist(x = Embeddings(seu)[idx, 1:30])
  for(j in (4:20)/20){
    tmp = silhouette(x = as.numeric(seu@meta.data[[
      paste0("integrated_snn_res.", j)]])[idx], dist = tmpD)
    tmp = data.table(tissue = i, resolution = j, 
                     silScore = mean(summary(tmp)$clus.avg.widths))
    oupSil = rbindlist(list(oupSil, tmp))
  }
  tmpD = oupSil[tissue == i & resolution >= 0.3][which.max(silScore)]
  p1 = ggplot(oupSil[tissue == i], 
              aes(resolution, silScore, label = round(silScore, 3))) + 
    geom_point() + geom_text(vjust = 1.35, size = 4) + 
    geom_point(data = tmpD, shape = 1, color = "red", size = 6) + 
    ylab("silhouette score") + theme_classic(base_size = 20)
  ggsave(p1, height = 6, width = 6, 
         filename = paste0("res_1inteData/mergeClust/", i, "_silScore.png"))
  
  # Generate pseudobulk and cor
  seu$seurat_clusters <- seu@meta.data[[
    paste0("integrated_snn_res.", tmpD$resolution)]]
  seu$seurat_clusters = factor(seu$seurat_clusters, levels = 
                                 sort(as.numeric(levels(seu$seurat_clusters))))
  tmp = matrix(nrow = uniqueN(seu$seurat_clusters), ncol = 30)
  rownames(tmp) = levels(seu$seurat_clusters)
  for(j in rownames(tmp)){
    tmp[j, ] = colMeans(Embeddings(seu)[seu$seurat_clusters == j, 1:30])
  }
  tmpCor <- cor(t(tmp), method = 'pearson') 
  tmpSize = 18; if(uniqueN(seu$seurat_clusters) > 15){tmpSize = 10}
  pheatmap::pheatmap(tmpCor, breaks = hmapBrk, color = hmapCol, border_color = "white", 
                     display_numbers = T, number_format = "%.02f", number_color = "black",
                     fontsize_row = tmpSize, fontsize_col = tmpSize,
                     fontsize_number = tmpSize, height = 8, width = 8, angle_col = 45, 
                     filename = paste0("res_1inteData/mergeClust/", i, "_spearman.png"))
  
  # Merge clusters and rename
  seu$merged_clusters = as.character(seu$seurat_clusters)
  if(i == "liver"){
    seu@meta.data[seu$seurat_clusters == "3", "merged_clusters"] = "3_5"
    seu@meta.data[seu$seurat_clusters == "5", "merged_clusters"] = "3_5"}
  if(i == "lung"){
    seu@meta.data[seu$seurat_clusters == "3", "merged_clusters"] = "3_12"
    seu@meta.data[seu$seurat_clusters =="12", "merged_clusters"] = "3_12"
    seu@meta.data[seu$seurat_clusters == "4", "merged_clusters"] = "4_8_13"
    seu@meta.data[seu$seurat_clusters == "8", "merged_clusters"] = "4_8_13"
    seu@meta.data[seu$seurat_clusters =="13", "merged_clusters"] = "4_8_13"
    seu@meta.data[seu$seurat_clusters == "1", "merged_clusters"] = "1_15"
    seu@meta.data[seu$seurat_clusters =="15", "merged_clusters"] = "1_15"}
  if(i == "heart"){
    seu@meta.data[seu$seurat_clusters == "4", "merged_clusters"] = "4_10"
    seu@meta.data[seu$seurat_clusters =="10", "merged_clusters"] = "4_10"}
  if(i == "endo"){
    seu@meta.data[seu$seurat_clusters == "2", "merged_clusters"] = "2_9"
    seu@meta.data[seu$seurat_clusters == "9", "merged_clusters"] = "2_9"
    seu@meta.data[seu$seurat_clusters == "1", "merged_clusters"] = "1_4"
    seu@meta.data[seu$seurat_clusters == "4", "merged_clusters"] = "1_4"
    seu@meta.data[seu$seurat_clusters == "0", "merged_clusters"] = "0_14"
    seu@meta.data[seu$seurat_clusters =="14", "merged_clusters"] = "0_14"
    seu@meta.data[seu$seurat_clusters == "8", "merged_clusters"] = "8_16"
    seu@meta.data[seu$seurat_clusters =="16", "merged_clusters"] = "8_16"}
  if(i == "kidney"){
    seu@meta.data[seu$seurat_clusters == "1", "merged_clusters"] = "1_6"
    seu@meta.data[seu$seurat_clusters == "6", "merged_clusters"] = "1_6"}
  seu$merged_clusters = factor(seu$merged_clusters)
  
  # Find marker genes for merged clusters
  Idents(seu) <- seu$merged_clusters
  tmpMarkers = FindAllMarkers(seu, only.pos = TRUE, 
                              logfc.threshold = 1, min.pct = 0.25)
  tmpMarkers = data.table(tmpMarkers)
  spp1Clust = tmpMarkers[gene == "SPP1"]$cluster
  rna1Clust = tmpMarkers[gene == "RNASE1"]$cluster
  fcn1Clust = tmpMarkers[gene == "FCN1"]$cluster
  
  # Annotate SPP1 cluster
  seu$clust = as.character(seu$merged_clusters)
  seu@meta.data[seu$merged_clusters %in% fcn1Clust, "clust"] = paste0(
    "fcn1-", seu@meta.data[seu$merged_clusters %in% fcn1Clust, "clust"])
  seu@meta.data[seu$merged_clusters %in% rna1Clust, "clust"] = paste0(
    "rnase-", seu@meta.data[seu$merged_clusters %in% rna1Clust, "clust"])
  seu@meta.data[seu$merged_clusters %in% spp1Clust, "clust"] = paste0(
    "spp1-", seu@meta.data[seu$merged_clusters %in% spp1Clust, "clust"])
  seu$clust = factor(seu$clust)
  Idents(seu) <- seu$clust
  
  # Add to seu
  tmpMarkers$clust = tmpMarkers$cluster
  tmpMarkers[cluster %in% fcn1Clust]$clust = paste0(
    "fcn1-", tmpMarkers[cluster %in% fcn1Clust]$clust)
  tmpMarkers[cluster %in% rna1Clust]$clust = paste0(
    "rnase-", tmpMarkers[cluster %in% rna1Clust]$clust)
  tmpMarkers[cluster %in% spp1Clust]$clust = paste0(
    "spp1-", tmpMarkers[cluster %in% spp1Clust]$clust)
  tmpMarkers = tmpMarkers[order(clust)]
  seu@misc$FindAllMarkers = tmpMarkers
  seu@misc$spp1Clust = spp1Clust
  seu@misc$rna1Clust = rna1Clust
  seu@misc$fcn1Clust = fcn1Clust
  print(tmpMarkers[gene == "SPP1"])
  print(tmpMarkers[gene == "RNASE1"])
  print(tmpMarkers[gene == "FCN1"])
  
  # Plot heatmap
  tmpC = c(); tmpCbrk = c()
  for(j in unique(seu$clust)){
    set.seed(42)
    tmptmp = colnames(seu)[seu$clust == j]
    tmpC = c(tmpC, sample(tmptmp, min(length(tmptmp), 200)))
    tmpCbrk = c(tmpCbrk, rep(j, min(length(tmptmp), 200)))
  }
  tmpG = unique(tmpMarkers[avg_log2FC > 3 | gene == "SPP1"]$gene)
  tmpG = tmpG[tmpG %in% rownames(seu@assays$integrated@scale.data)]
  ggData = seu@assays$integrated@scale.data[tmpG, tmpC]
  ggData[ggData >  2] =  2; ggData[ggData < -2] = -2
  ha = rowAnnotation(gene = anno_mark(
    at = which(tmpG %in% genes.to.label), 
    labels = tmpG[which(tmpG %in% genes.to.label)]))
  hb = HeatmapAnnotation(clust = tmpCbrk, annotation_name_side = "left")
  png(file = paste0("res_1inteData/mergeClust/", i, "_markers.png"), 
      width = 8, height = 8, res = 300, units = "in")
  draw(Heatmap(ggData, name = " ", column_split = tmpCbrk,
               top_annotation = hb, right_annotation = ha, 
               cluster_rows = FALSE, cluster_columns = FALSE, 
               show_column_names = FALSE, show_row_names = FALSE))
  dev.off()
  
  # Extra plot
  p1 = DimPlot(seu, group.by = "clust", label = TRUE, raster = FALSE,
               cols = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$clust)))
  p2 = DotPlot(seu, group.by = "clust", 
               features = c("RNASE1", "FCN1", "SPP1", "TREM2", "MMP9")) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds")))
  ggsave(p1 + p2, height = 6, width = 14, 
         filename = paste0("res_1inteData/mergeClust/", i, "_clustGEX.png"))
  
  
  ### 4. Final Plot
  # Plot 
  colIdent = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$orig.ident))
  ggOut = list()
  ggOut[["sample"]] = DimPlot(seu, group.by = "orig.ident", cols = colIdent, 
                              shuffle = T, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  ggOut[["clust"]] = DimPlot(seu, group.by = "clust", label = T, cols = 
                               colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$clust)), 
                              shuffle = T, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  tmpG = c("SPP1", "MMP9", "TREM2", "RNASE1")
  for(j in tmpG){
    if(j %in% rownames(seu)){DefaultAssay(seu) = "integrated"}
    if(!j %in% rownames(seu)){DefaultAssay(seu) = "RNA"}
    ggOut[[j]] = FeaturePlot(seu, features = j, order = T, pt.size = 0.25) + 
      scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  ggsave(wrap_plots(ggOut), height = 12, width = 16,
         filename = paste0("res_1inteData/", i, "_umap.png"))
  saveRDS(seu, file = paste0("res_1inteData/", i, ".rds"))
}
fwrite(oupSil, sep = "\t", file = "res_1inteData/silScore.txt")




