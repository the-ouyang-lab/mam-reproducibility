rm(list = ls())
library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(cluster)
library(ComplexHeatmap)
library(msigdbr)
library(clusterProfiler)
if(!dir.exists("res_1inteData/")){dir.create("res_1inteData/")}
if(!dir.exists("res_1inteData/comet/")){dir.create("res_1inteData/comet/")}

hmapCol <- colorRampPalette(c("blue", "white", "red"))(21)
hmapBrk <- seq(-1, 1, length.out = 21)
genes.to.label <- c(
  'SELENOP','C1QA','RNASE1' ,'C1QC','FCN1', 'SDC2', 'PLA2G7','S100A9','S100A8',
  'VCAN', 'APOE','GPNMB','CD9','STAB1','SPP1', 'FABP5', 'TREM2', 'F13A1', 
  'CCL20', 'FOLR2', 'INHBA',  'CHI3L1',  'MATK','PTGS2','IL1B', 'FABP4', 'LYVE1', 
  'APOC1', 'CCL4','S100A12', 'CTSK', 'LST1', 'CDKN1C', 'FCGR3A',
  'CCL3', 'DAB2','CCL4L2', 'PHLDA3', 'EGR1', 'ATF3', 'VEGFA', 'SEMA6B', 
  'IFI44L', 'IFIT3', 'IFIT1', 'MT1G', 'MT1X', 'MT1H')
samGene <- c("CD9", "TREM2", "SPP1")
mamGene <- fread("res_1inteData/mam_DEsumm.txt")

# liver   10401
# lung   184297
# heart   24461
# skin     2776
# endo     6030
# kidney   7965



### 1. Start Integration with TAM/LAM/SAM
## A. IO each dataset
# Note: SAMs are defined internally from liver
seuSPP1 <- readRDS("res_1inteData/SPP1seu.rds") 
seuTAM <- readRDS(paste0("~/../../publicData/macroVerse/", 
                         "2021_MoMac_VERSE.RDS"))  # TAM is cluster 3
seuLAM <- readRDS("rawData/0mac/myel_lams.RDS")    # LAM is cluster 5
seuDAM <- readRDS("rawData/0mac/dams_object.RDS")  # DAM is cluster 1/4/12
seuLiv <- readRDS("res_1inteData/liver.rds")       # to derive SAMs


## B. Derive SAMs
seuLiv <- AddModuleScore(seuLiv, features = list(sam = samGene), name = "samSig")
seuLiv$samCells = FALSE
seuLiv@meta.data[scale(seuLiv$samSig1) > 1.96, "samCells"] = TRUE
seuLiv$mamCells = FALSE
seuLiv@meta.data[intersect(colnames(seuSPP1)[seuSPP1$mam == "MAM"], colnames(seuLiv)), "mamCells"] = TRUE
p1 = DimPlot(seuLiv, group.by = "samCells", label = T, pt.size = 1, cols = c("grey","firebrick"))
p2 = DimPlot(seuLiv, group.by = "mamCells", label = T, pt.size = 1, cols = c("grey","firebrick"))
ggsave(p1 + p2, height = 6, width = 16, filename = paste0("res_1inteData/SAMs_umap.png"))


## C. Generate files for cometSC
tmp = as.matrix(seuSPP1@assays$SCT@data)
write.table(tmp, sep = "\t", quote = FALSE,
            file = "res_1inteData/comet/tabmarker.txt")
tmp = seuSPP1@reductions$umap@cell.embeddings
write.table(tmp, sep = "\t", quote = FALSE, col.names = FALSE,
            file = "res_1inteData/comet/tabvis.txt")
tmp = data.frame(clust = seuSPP1$mam, row.names = colnames(seuSPP1))
write.table(tmp, sep = "\t", quote = FALSE, col.names = FALSE,
            file = "res_1inteData/comet/tabcluster.txt")



### 2. Start integration
## A. prepare GEX
inpGEX = list()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  tmp = subset(seuSPP1, tissue == i)
  inpGEX[[i]] = tmp@assays$RNA@counts
}

tmp = subset(seuTAM, Clusters == "TREM2 Macrophage -3")
inpGEX[["tam"]] = tmp@assays$RNA@counts[intersect(rownames(tmp@assays$RNA@counts), rownames(seuSPP1)), ]
colnames(inpGEX[["tam"]]) = paste0(tmp$Study, "_", colnames(inpGEX[["tam"]]))

tmp = subset(seuDAM, seurat_clusters %in% c(1,4,12))
inpGEX[["dam"]] = tmp@assays$RNA@counts[intersect(rownames(tmp), rownames(seuSPP1)), ]
colnames(inpGEX[["dam"]]) = paste0(tmp$sampleName, "_", colnames(inpGEX[["dam"]]))

tmp = subset(seuLAM, seurat_clusters %in% c(5))
inpGEX[["lam"]] = tmp@assays$RNA@counts[intersect(rownames(tmp), rownames(seuSPP1)), ]
colnames(inpGEX[["lam"]]) = paste0(tmp$Patient, "_", colnames(inpGEX[["lam"]]))

inpGEX[["sam"]] = seuLiv@assays$RNA@counts[intersect(rownames(
  seuLiv@assays$RNA@counts), rownames(seuSPP1)), seuLiv$samCells == TRUE]
colnames(inpGEX[["sam"]]) = paste0("sam", colnames(inpGEX[["sam"]]))

# Start SCTransform
oupSeuList = list()
for(i in names(inpGEX)){
  tmp = inpGEX[[i]]
  tmp = CreateSeuratObject(counts = tmp, project = i)
  tmp = SCTransform(tmp, vars.to.regress = "orig.ident", verbose = F)
  tmp$tissue = i
  oupSeuList[[i]] = tmp
}

# Start integration
features <- SelectIntegrationFeatures(object.list = oupSeuList, nfeatures = 3000)
seu <- PrepSCTIntegration(object.list = oupSeuList, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seu, normalization.method = "SCT",
                                  anchor.features = features)
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT", k.weight = 50)
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, verbose = FALSE)

# Add some annotations
seu$group = "nil"
seu@meta.data[colnames(seuSPP1)[seuSPP1$mam == "MAM"], "group"]    = "SPP1+MAM+"
seu@meta.data[colnames(seuSPP1)[seuSPP1$mam == "nonMAM"], "group"] = "SPP1+MAM-"
seu@meta.data[seu$tissue == "tam", "group"] = "TAM"
seu@meta.data[seu$tissue == "lam", "group"] = "LAM"
seu@meta.data[seu$tissue == "dam", "group"] = "DAM"
seu@meta.data[seu$tissue == "sam", "group"] = "SAM"
seu$group = factor(seu$group, levels = c(
  "SPP1+MAM+","SPP1+MAM-","TAM","LAM","DAM","SAM"))
p1 = DimPlot(seu, group.by = "group", cols = c(
  "firebrick","darkorange","dodgerblue","purple","olivedrab","steelblue"))
ggsave(p1, height = 6, width = 8, filename = paste0("res_1inteData/allMac_umap.png"))

# Calculate scores
DefaultAssay(seu) <- "SCT"
seu <- AddModuleScore(seu, features = list(mamGene[mamGenes == TRUE]$gene))
seu$mamScore = seu$Cluster1; seu$Cluster1 = NULL
seu <- AddModuleScore(seu, features = list(mamGene[mamGenesECM == TRUE]$gene))
seu$mamECM = seu$Cluster1; seu$Cluster1 = NULL
seu <- AddModuleScore(seu, features = list(mamGene[mamGenesMetab == TRUE]$gene))
seu$mamMetab = seu$Cluster1; seu$Cluster1 = NULL

p1 <- VlnPlot(seu, features = c("mamECM", "mamMetab"), group.by = "group", pt.size = 0.001) &
  stat_summary(fun = median, geom = "point", size = 15, colour = "white", 
               shape = 95) & NoLegend()
ggsave(p1, height = 6, width = 10, filename = paste0("res_1inteData/allMac_score.png"))


### Final IO
saveRDS(seu, file = paste0("res_1inteData/AllMACseu.rds"))




