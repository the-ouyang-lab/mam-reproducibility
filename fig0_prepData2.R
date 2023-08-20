rm(list = ls())
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(Matrix)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(harmony)
library(msigdbr)
if(!dir.exists("res_0prepData2/")){dir.create("res_0prepData2/")}
if(!dir.exists("res_0prepData2rds/")){dir.create("res_0prepData2rds/")}



### Find common genes for each tissue then subset each matrix...
# round2 is exclusion
inpParams = data.table(
  study = c("liver1","liver2","lung1","lung2","lung3","lung4","heart1","heart2",
            "skin1","skin2","endo1","endo2","kidney1","kidney2","kidney3"),
  round2 = c("10|7|12|3|13|14", "11|8|2|0|12",                      # liver1/2
             "15|9|16|12|17|6", "18|7|16|11|10|17",                 # lung1/2
             "18|13|9", "13|9|16|10",                               # lung3/4
             "14|6|13|11|9", "15|12|8|3",                           # heart1/2
             "1|7", "10|0|4|3|2|9|6|8|5",                           # skin1/2
             "12|16|19|13|9|18|11|15|17", "16|8|7|2|14|1|5|11|13|12|8|15|9", # endo1/2
             "15|11|16|7|5|13|14|9|8","13|10|7|3|9|15|6|12|11|8",   #kidney1-3
             "10|6|3|12|9|8|11") 
)
filterGenes = c("SPP1", "MMP7", "MMP9",
                "GZMB", "GNLY", "CCR7",     # NK
                "CD1C", "FCER1A",           # cDC
                "STMN1", "TUBB",            # prolif
                "DCN", "LUM")               # mesenchymal
# Can also use SFTPB for lung3
tmp42 = data.table()
for(j in inpParams$study){
  tmp = fread(paste0("res_0prepData1/", j, "_round1.txt"))
  tmp = tmp[features.plot == "CD68"][avg.exp.scaled > 0][pct.exp > 10]
  tmp$dataset = j  
  tmp42 = rbindlist(list(tmp42, tmp))
}



### Count num cells and genes in each dataset
oupNum = data.table()
for(j in c("a","b")){
  tmp = fread(paste0("res_0prepData1/oupNUMr2", j, ".txt"))
  oupNum = rbindlist(list(oupNum, tmp))
}
oupNum$r3nCell = 42
oupNum$r3nGene = 42
gc()


# Start cleaning 
for(j in inpParams$study){
  inpGEX = readRDS(paste0("res_0prepData1rds/oupGEXr2_", j, ".rds"))
  seu = CreateSeuratObject(counts = inpGEX, project = j)
  seu = NormalizeData(seu, verbose=F)
  seu = FindVariableFeatures(seu, verbose=F)
  seu = ScaleData(seu, verbose=F)
  seu = RunPCA(seu, verbose=F)
  if(j %in% c("skin1", "skin2")){
    seu = RunUMAP(seu, dims = 1:30, verbose=F)
    seu = FindNeighbors(seu, dims = 1:30, verbose=F)
  } else {
    seu = RunHarmony(seu, group.by.vars = "orig.ident", project.dim=F, verbose=F)
    seu = RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose=F)
    seu = FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose=F)
  }
  seu = FindClusters(seu, resolution = 1.0, verbose=F)
  seu$pctMT = 100 * colSums(inpGEX[grep("^MT-", rownames(inpGEX)),]) / colSums(inpGEX)
  
  colIdent = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$orig.ident))
  colClust = colorRampPalette(brewer.pal(8, "Dark2"))(uniqueN(seu$seurat_clusters))
  ggOut = list()
  ggOut[["sample"]] = DimPlot(seu, group.by = "orig.ident", cols = colIdent, 
                              shuffle = T, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  ggOut[["clust"]] = DimPlot(seu, group.by = "seurat_clusters", label = T, cols = colClust, 
                             label.size = 3, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  names(colClust) = levels(Idents(seu))
  colClust[strsplit(inpParams[study == j]$round2, "\\|")[[1]]] = "grey"
  ggOut[["clust2"]] = DimPlot(seu, group.by = "seurat_clusters", label = T, cols = colClust, 
                             label.size = 3, pt.size = 0.25) + NoLegend() + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
  for(i in filterGenes){
    if(!i %in% rownames(seu)){i = "CD68"}
    ggOut[[i]] = FeaturePlot(seu, features = i, order = T, pt.size = 0.25) + 
      scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed() +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }
  ggOut[["dotPlot"]] = DotPlot(seu, features = filterGenes[filterGenes %in% rownames(seu)]) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + 
    theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust=0))
  layout <- "ABCPPPP
             DEFPPPP
             GHIPPPP
             JKLPPPP
             MNOPPPP"
  ggsave(wrap_plots(ggOut) + plot_layout(design = layout), height = 10, width = 15,
         filename = paste0("res_0prepData2/", j, "_round2.png"))
  fwrite(data.table(ggOut[["dotPlot"]]$data)[order(-features.plot, -pct.exp)], sep = "\t",
         file = paste0("res_0prepData2/", j, "_round2.txt"))

  ## Subset
  tmp = strsplit(inpParams[study == j]$round2, "\\|")[[1]]
  seu = subset(seu, idents = tmp, invert = TRUE)
  inpGEX = inpGEX[, colnames(seu)]
  inpGEX = inpGEX[rowSums(inpGEX) > 0, ]
  oupNum[study == j]$r3nCell = ncol(inpGEX)
  oupNum[study == j]$r3nGene = nrow(inpGEX)
  fwrite(oupNum, sep = "\t", file = "res_0prepData2/oupNUMr3b.txt")
  saveRDS(inpGEX, file = paste0("res_0prepData2rds/oupGEXr3_", j, ".rds"))
  gc()
}




