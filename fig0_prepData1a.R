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
if(!dir.exists("res_0prepData1/")){dir.create("res_0prepData1/")}
if(!dir.exists("res_0prepData1rds/")){dir.create("res_0prepData1rds/")}
oupGEX = list()


### Gather and compile counts matrix
## Liver_Ramachandran_GSE136103
inpFiles = list.files("rawData/Liver_Ramachandran_GSE136103/", pattern = ".mtx.gz")
inpFiles = gsub("matrix.mtx.gz", "", inpFiles)
oupGEX$liver1 = NULL
for(i in inpFiles){
  tmp = ReadMtx(mtx = paste0("rawData/Liver_Ramachandran_GSE136103/", i, "matrix.mtx.gz"),
                cells = paste0("rawData/Liver_Ramachandran_GSE136103/", i, "barcodes.tsv.gz"),
                features = paste0("rawData/Liver_Ramachandran_GSE136103/", i, "genes.tsv.gz"))
  tmpID = paste0("liver1-", strsplit(i, "_")[[1]][2], "_", strsplit(i, "_")[[1]][3])
  tmpID = gsub("healthy", "norm-", tmpID)
  tmpID = gsub("cirrhotic", "cirr-", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$liver1)){
    oupGEX$liver1 = tmp
  } else {
    oupGEX$liver1 = cbind(oupGEX$liver1, tmp)
  }
}


## Liver_Fred_authorData
oupGEX$liver2 = NULL
for(i in 1:10){
  tmp = ReadMtx(mtx = paste0("rawData/Liver_Fred_authorData/", i, "/matrix.mtx"),
                cells = paste0("rawData/Liver_Fred_authorData/", i, "/barcodes.tsv"),
                features = paste0("rawData/Liver_Fred_authorData/", i, "/genes.tsv"))
  tmpID = paste0("liver2-NASH-patient", i)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$liver2)){
    oupGEX$liver2 = tmp
  } else {
    commonGenes = intersect(rownames(oupGEX$liver2), rownames(tmp))
    oupGEX$liver2 = cbind(oupGEX$liver2[commonGenes, ], tmp[commonGenes, ])
  }
}
gc()


## Lung_Morse_GSE128033
# https://erj.ersjournals.com/content/54/2/1802441.figures-only#fig-data-supplementary-materials
# Here, we exclude SC14, SC31, SC31D due to V1 and low cell number
# Here, we exclude SC228 and SC249 due to bal sample
# Here, we exclude SC89 and SC95 due to mac depletion
inpFiles = list.files("rawData/Lung_Morse_GSE128033/", pattern = ".mtx.gz")
inpFiles = inpFiles[c(4:8,11,12,14,15,17,18)]
inpFiles = gsub("matrix.mtx.gz", "", inpFiles)
oupGEX$lung1 = NULL
for(i in inpFiles){
  tmp = ReadMtx(mtx = paste0("rawData/Lung_Morse_GSE128033/", i, "matrix.mtx.gz"),
                cells = paste0("rawData/Lung_Morse_GSE128033/", i, "barcodes.tsv.gz"),
                features = paste0("rawData/Lung_Morse_GSE128033/", i, "genes.tsv.gz"))
  tmp = tmp[, colSums(tmp) > 300]
  tmpID = strsplit(i, "_")[[1]][2]
  if(grepl("IPF", tmpID)){
    tmpID = paste0("lung1-IPF-", gsub("IPF|UP|LOW","", tmpID))
  } else {
    tmpID = paste0("lung1-norm-", gsub("NOR|UP|LOW","", tmpID))
  }
  tmpID = gsub("SC155", "SC155_LOW", tmpID)
  tmpID = gsub("SC156", "SC155_UPP", tmpID)
  tmpID = gsub("SC87",  "SC87_LOW", tmpID)
  tmpID = gsub("SC88",  "SC87_UPP", tmpID)
  tmpID = gsub("SC93",  "SC93_LOW", tmpID)
  tmpID = gsub("SC94",  "SC93_UPP", tmpID)
  tmpID = gsub("SC153", "SC153_LOW", tmpID)
  tmpID = gsub("SC154", "SC153_UPP", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$lung1)){
    oupGEX$lung1 = tmp
  } else {
    oupGEX$lung1 = cbind(oupGEX$lung1, tmp)
  }
}


## Lung_Reyfman_GSE122960
# Here, we remove cryobiopsy
inpFiles = list.files("rawData/Lung_Reyfman_GSE122960/", pattern = ".h5")[-5]
oupGEX$lung2 = NULL
for(i in inpFiles){
  tmp = Read10X_h5(paste0("rawData/Lung_Reyfman_GSE122960/", i))
  tmpID = paste0("lung2-", strsplit(i, "_")[[1]][2], "-", 
                 gsub("0", "", strsplit(i, "_")[[1]][3]))
  tmpID = gsub("Donor", "norm", tmpID)
  tmpID = gsub("Myositis-", "myositis", tmpID)
  tmpID = gsub("SSc-", "SSc", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$lung2)){
    oupGEX$lung2 = tmp
  } else {
    oupGEX$lung2 = cbind(oupGEX$lung2, tmp)
  }
}

## Lung_Adams_GSE136831
tmp = ReadMtx(mtx = "rawData/Lung_Adams_GSE136831/GSE136831_RawCounts_Sparse.mtx.gz",
              cells = "rawData/Lung_Adams_GSE136831/GSE136831_AllCells.cellBarcodes.txt.gz",
              features = "rawData/Lung_Adams_GSE136831/GSE136831_AllCells.GeneIDs.txt.gz", 
              skip.feature = 1)
tmpMeta = fread("rawData/Lung_Adams_GSE136831/GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
# all.equal(colnames(tmp), tmpMeta$CellBarcode_Identity)
tmpMeta = tmpMeta[Disease_Identity != "COPD"]
tmpMeta = tmpMeta[Manuscript_Identity %in% c('Macrophage_Alveolar', 'Macrophage', 'cMonocyte')]
tmpMeta[Disease_Identity == "Control"]$Disease_Identity = "norm"
tmpMeta = tmpMeta[!Subject_Identity %in% c("001C","048C","098C","133C",
                                           "1372C","137C","160C","192C",
                                           "222C","244C","388C","454C")] # CHECK
tmpMeta$tmpID = paste0("lung3-", tmpMeta$Disease_Identity, "-", tmpMeta$Subject_Identity, "_",
                       tmpMeta$CellBarcode_Identity)
# sort(unique(paste0("lung3-", tmpMeta$Disease_Identity, "-", tmpMeta$Subject_Identity)))
oupGEX$lung3 = tmp[, tmpMeta$CellBarcode_Identity]
colnames(oupGEX$lung3) = tmpMeta$tmpID
for(i in unique(tstrsplit(colnames(oupGEX$lung3), "_")[[1]])){
  print(paste0(i, " | ", dim(oupGEX$lung3[, grep(i, colnames(oupGEX$lung3))])[1], 
                  " | ", dim(oupGEX$lung3[, grep(i, colnames(oupGEX$lung3))])[2]))
}


## Lung_Valenzi_GSE128169_GSE156310
# https://www.frontiersin.org/articles/10.3389/fimmu.2021.595811/full
# Here we use GSE128169_RAW.tar: SC51,SC52,SC63,SC64,SC108,SC109,SC135,SC136
# Here we use GSE156310_RAW.tar: SC174, SC175
inpFiles = list.files("rawData/Lung_Valenzi_GSE128169_GSE156310/", pattern = ".mtx.gz")
inpFiles = inpFiles[6:13]
inpFiles = gsub("matrix.mtx.gz", "", inpFiles)
oupGEX$lung4 = NULL
for(i in inpFiles){
  tmp = ReadMtx(mtx = paste0("rawData/Lung_Valenzi_GSE128169_GSE156310/", i, "matrix.mtx.gz"),
                cells = paste0("rawData/Lung_Valenzi_GSE128169_GSE156310/", i, "barcodes.tsv.gz"),
                features = paste0("rawData/Lung_Valenzi_GSE128169_GSE156310/", i, "genes.tsv.gz"))
  tmp = tmp[, colSums(tmp) > 300]
  tmpID = strsplit(i, "_")[[1]][2]
  if(grepl("SSC", tmpID)){
    tmpID = paste0("lung4-SSC-", gsub("SSC|UP|LOW","", tmpID))
  } else {
    tmpID = paste0("lung4-norm-", gsub("NOR|UP|LOW","", tmpID))
  }
  tmpID = gsub("SC51",  "SC51_LOW", tmpID)
  tmpID = gsub("SC52",  "SC51_UPP", tmpID)
  tmpID = gsub("SC63",  "SC63_LOW", tmpID)
  tmpID = gsub("SC64",  "SC63_UPP", tmpID)
  tmpID = gsub("SC108", "SC108_LOW", tmpID)
  tmpID = gsub("SC109", "SC108_UPP", tmpID)
  tmpID = gsub("SC135", "SC135_LOW", tmpID)
  tmpID = gsub("SC136", "SC135_UPP", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$lung4)){
    oupGEX$lung4 = tmp
  } else {
    oupGEX$lung4 = cbind(oupGEX$lung4, tmp)
  }
}
tmp = Read10X_h5("rawData/Lung_Valenzi_GSE128169_GSE156310/GSM4728858_SC174raw_feature_bc_matrix.h5")
tmp = tmp[, colSums(tmp) > 300]
tmpID = "lung4-IPF-SC174_LOW"; print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
commonGenes = intersect(rownames(oupGEX$lung4), rownames(tmp))
oupGEX$lung4 = cbind(oupGEX$lung4[commonGenes, ], tmp[commonGenes, ])
tmp = Read10X_h5("rawData/Lung_Valenzi_GSE128169_GSE156310/GSM4728859_SC175raw_feature_bc_matrix.h5")
tmp = tmp[, colSums(tmp) > 300]
tmpID = "lung4-IPF-SC174_UPP"; print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
commonGenes = intersect(rownames(oupGEX$lung4), rownames(tmp))
oupGEX$lung4 = cbind(oupGEX$lung4[commonGenes, ], tmp[commonGenes, ])
gc()


## Heart_Koenig_GSE183852
load("rawData/Heart_Koenig_GSE183852/GSE183852_DCM_Cells.Robj")
oupGEX$heart1 = HDCM@assays$RNA@counts
colnames(oupGEX$heart1) = gsub("HDCM", "heart1-DCM-HDCM", colnames(oupGEX$heart1))
colnames(oupGEX$heart1) = gsub("heart1-DCM-HDCM5", "heart1-norm-HDCM5", colnames(oupGEX$heart1))
colnames(oupGEX$heart1) = gsub("heart1-DCM-HDCM7", "heart1-norm-HDCM7", colnames(oupGEX$heart1))
for(i in unique(tstrsplit(colnames(oupGEX$heart1), "_")[[1]])){
  print(paste0(i, " | ", dim(oupGEX$heart1[, grep(i, colnames(oupGEX$heart1))])[1], 
                  " | ", dim(oupGEX$heart1[, grep(i, colnames(oupGEX$heart1))])[2]))
}
rm(HDCM)


## Heart_Rao_GSE145154
inpFiles = list.files("rawData/Heart_Rao_GSE145154/", pattern = ".mtx.gz")
inpFiles = inpFiles[!grepl("Bld|N_|_RM", inpFiles)]
inpFiles = gsub("matrix.mtx.gz", "", inpFiles)
oupGEX$heart2 = NULL
for(i in inpFiles){
  tmp = ReadMtx(mtx = paste0("rawData/Heart_Rao_GSE145154/", i, "matrix.mtx.gz"),
                cells = paste0("rawData/Heart_Rao_GSE145154/", i, "barcodes.tsv.gz"),
                features = paste0("rawData/Heart_Rao_GSE145154/", i, "features.tsv.gz"))
  tmpID = paste0("heart2-", strsplit(i, "_")[[1]][2])
  tmpID = gsub("-N-", "-norm-", tmpID)
  tmpID = gsub("-L", "_L", tmpID)
  tmpID = gsub("-R", "_R", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$heart2)){
    oupGEX$heart2 = tmp
  } else {
    oupGEX$heart2 = cbind(oupGEX$heart2, tmp)
  }
}
gc()



### Count num cells and genes in each dataset
oupNum = data.table()
for(j in names(oupGEX)){
  oupGEX[[j]] = oupGEX[[j]][rowSums(oupGEX[[j]]) > 0, ]
  tmp = data.table(study = j,
                   r0nCell = ncol(oupGEX[[j]]),
                   r0nGene = nrow(oupGEX[[j]]))
  oupNum = rbindlist(list(oupNum, tmp))
}
oupNum$r1nCell = 42
oupNum$r1nGene = 42
oupNum$r2nCell = 42
oupNum$r2nGene = 42
gc()


### Subset Myeloid cells
for(j in names(oupGEX)){
  # Start cellQC for each dataset
  oupMeta = data.table(cell = colnames(oupGEX[[j]]))
  oupMeta$nGene = colSums(oupGEX[[j]] > 0)
  oupMeta$pctMT = 100 * colSums(oupGEX[[j]][
    grep("^MT-", rownames(oupGEX[[j]])),]) / colSums(oupGEX[[j]])
  oupMeta$pctHB = 100 * colSums(oupGEX[[j]][
    grep("^HBB|^HBA", rownames(oupGEX[[j]])),]) / colSums(oupGEX[[j]])
  oupMeta = oupMeta[nGene > 300 & nGene < 5000]
  oupMeta = oupMeta[pctMT < 20 & pctHB < 0.1]
  oupGEX[[j]] = oupGEX[[j]][, oupMeta$cell]
  oupGEX[[j]] = oupGEX[[j]][rowSums(oupGEX[[j]]) > 0, ]
  oupNum[study == j]$r1nCell = ncol(oupGEX[[j]])
  oupNum[study == j]$r1nGene = nrow(oupGEX[[j]])
  
  # Start first round Seurat to find CD45+ CD68+
  seu = CreateSeuratObject(counts = oupGEX[[j]], project = j)
  seu = NormalizeData(seu, verbose=F)
  seu = FindVariableFeatures(seu, verbose=F)
  seu = ScaleData(seu, verbose=F)
  seu = RunPCA(seu, verbose=F)
  seu = RunHarmony(seu, group.by.vars = "orig.ident", project.dim=F, verbose=F)
  seu = RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose=F)
  seu = FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose=F)
  seu = FindClusters(seu, resolution = 1.0, verbose=F)
  
  colIdent = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$orig.ident))
  colClust = colorRampPalette(brewer.pal(8, "Dark2"))(uniqueN(seu$seurat_clusters))
  p1 = DimPlot(seu, group.by = "orig.ident", cols = colIdent, shuffle = T) + 
    NoLegend() + coord_fixed()
  p2 = DimPlot(seu, group.by = "seurat_clusters", label = T, cols = colClust, 
               label.size = 5) + NoLegend() + coord_fixed()
  p3 = FeaturePlot(seu, features = c("PTPRC"), order = T) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed()
  tmp4 = "CD68"; if(j == "lung3"){tmp4 = "CD14"}
  p4 = FeaturePlot(seu, features = c(tmp4), order = T) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed()
  p5 = FeaturePlot(seu, features = c("SPP1"), order = T) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed()
  p6 = FeaturePlot(seu, features = c("MMP7"), order = T) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds"))) + coord_fixed()
  
  p7 = DotPlot(seu, features = c("PTPRC", tmp4, "SPP1", "MMP7")) + 
    scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds")))
  layout <- "ABG
             CDG
             EFG"
  ggsave(p1 + p2 + p3 + p4 + p5 + p6 + p7 + plot_layout(design = layout), 
         height = 12, width = 16,
         filename = paste0("res_0prepData1/", j, "_round1.png"))
  fwrite(data.table(p7$data)[order(-features.plot, -pct.exp)], sep = "\t",
         file = paste0("res_0prepData1/", j, "_round1.txt"))
  
  ## Subset
  if(j == "liver2"){
    tmp = data.table(p7$data)[features.plot == "CD68" & pct.exp > 25 &
                                avg.exp.scaled > 0]$id; print(tmp)
    seu = subset(seu, idents = tmp)
  } else if(j == "heart2"){
      tmp = data.table(p7$data)[features.plot == "CD68" & pct.exp > 30 &
                                  avg.exp.scaled > 0]$id; print(tmp)
      seu = subset(seu, idents = tmp)
  } else if(j != "lung3"){
    tmp = data.table(p7$data)[features.plot == "CD68" & pct.exp > 40 &
                                avg.exp.scaled > 0]$id; print(tmp)
    seu = subset(seu, idents = tmp)
  }
  oupGEX[[j]] = oupGEX[[j]][, colnames(seu)]
  oupGEX[[j]] = oupGEX[[j]][rowSums(oupGEX[[j]]) > 0, ]
  oupNum[study == j]$r2nCell = ncol(oupGEX[[j]])
  oupNum[study == j]$r2nGene = nrow(oupGEX[[j]])
  fwrite(oupNum, sep = "\t", file = "res_0prepData1/oupNUMr2a.txt")
  saveRDS(oupGEX[[j]], file = paste0("res_0prepData1rds/oupGEXr2_", j, ".rds"))
  gc()
}




