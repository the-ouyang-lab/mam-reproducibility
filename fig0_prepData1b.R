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
## Skin_Gur_GSE195452
tmpMeta = fread("rawData/Skin_Gur_GSE195452/GSE195452-GPL18573_series_matrix.txt",
                 skip = 32, fill = TRUE); tmpMeta = tmpMeta[1:11, ]
tmpMeta = data.table(sample = colnames(tmpMeta)[-1],
                     geo = as.character(tmpMeta[1, -1]),
                     source = as.character(tmpMeta[7, -1]),
                     tissue = gsub("^tissue: ", "", as.character(tmpMeta[9, -1])),
                     selection = gsub("^selection marker: ", "", as.character(tmpMeta[10, -1])),
                     patient = gsub("^patient id: ", "", as.character(tmpMeta[11, -1])))
tmpMeta$filename = paste0(tmpMeta$geo, "_", tmpMeta$sample, ".txt.gz")
tmpMeta$type = "Ctrl"; tmpMeta[grep("^pt", patient)]$type = "SSc"
tmpMeta = tmpMeta[tissue != "Blood"][selection == "CD45+"][patient != "-"]
# inpFiles = list.files("rawData/Skin_Gur_GSE195452/", pattern = ".txt.gz")
oupGEX$skin1 = NULL
for(i in 1:nrow(tmpMeta)){
  tmp = fread(paste0("rawData/Skin_Gur_GSE195452/", tmpMeta[i]$filename))
  tmp = as(as.matrix(data.frame(tmp[, -1], row.names = tmp$V1)), "dgCMatrix")
  tmpID = paste0("skin1-", tmpMeta[i]$patient)
  tmpID = gsub("-Ctrl", "-norm-Ctrl", tmpID)
  tmpID = gsub("-pt", "-SSc-pt", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$skin1)){
    oupGEX$skin1 = tmp
  } else {
    oupGEX$skin1 = cbind(oupGEX$skin1, tmp)
  }
}


## Skin_Deng_GSE163973
inpFiles = list.dirs("rawData/Skin_Deng_GSE163973/")[-1]
oupGEX$skin2 = NULL
for(i in inpFiles){
  tmp = Read10X(data.dir = i)
  tmpID = paste0("skin2-", strsplit(i, "//")[[1]][2])
  tmpID = gsub("_matrix", "", tmpID)
  tmpID = gsub("NF", "norm-NF", tmpID)
  tmpID = gsub("KF", "Keloid-KF", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$skin2)){
    oupGEX$skin2 = tmp
  } else {
    oupGEX$skin2 = cbind(oupGEX$skin2, tmp)
  }
}
gc()


## Endometrium_Tan_GSE179640
inpFiles = list.files("rawData/Endometrium_Tan_GSE179640/", pattern = ".h5")
inpFiles = inpFiles[grep("Ctrl|EuE", inpFiles)]
oupGEX$endo1 = NULL
for(i in inpFiles){
  tmp = Read10X_h5(paste0("rawData/Endometrium_Tan_GSE179640/", i))
  tmpID = paste0("endo1-", strsplit(i, "_")[[1]][3], "-",
                 strsplit(i, "_")[[1]][2], "_", strsplit(i, "_")[[1]][3])
  tmpID = gsub("Ctrl", "norm", tmpID)
  tmpID = gsub("EuE", "endo", tmpID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$endo1)){
    oupGEX$endo1 = tmp
  } else {
    commonGenes = intersect(rownames(oupGEX$endo1), rownames(tmp))
    oupGEX$endo1 = cbind(oupGEX$endo1[commonGenes, ], tmp[commonGenes, ])
  }
}


## Endometrium_Fonseca_GSE213216
# Require https://github.com/lawrenson-lab/AtlasEndometriosis/blob/main/GEO_id_patient_class.csv
tmpMeta = fread("rawData/Endometrium_Fonseca_GSE213216/GEO_id_patient_class.csv")
tmpMeta$type = "norm"
tmpMeta[`Major Class` %in% c("Endometrioma", "Endometriosis")]$type = "endo"
tmpMeta$ID = gsub("Sample", "sample", tmpMeta$ID)
tmpMeta = tmpMeta[-c(2,3,4,6), ]
tmpMeta = tmpMeta[`Major Class` %in% c("Endometriosis","Eutopic Endometrium")]
oupGEX$endo2 = NULL
for(i in 1:nrow(tmpMeta)){
  tmp = ReadMtx(mtx = paste0("rawData/Endometrium_Fonseca_GSE213216/", tmpMeta[i]$ID,
                             "/outs/filtered_feature_bc_matrix/matrix.mtx.gz"),
                cells = paste0("rawData/Endometrium_Fonseca_GSE213216/", tmpMeta[i]$ID,
                               "/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"),
                features = paste0("rawData/Endometrium_Fonseca_GSE213216/", tmpMeta[i]$ID,
                                  "/outs/filtered_feature_bc_matrix/features.tsv.gz"))
  tmpID = paste0("endo2-", tmpMeta[i]$type, "-pat", tmpMeta[i]$Patient, tmpMeta[i]$ID)
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$endo2)){
    oupGEX$endo2 = tmp
  } else {
    commonGenes = intersect(rownames(oupGEX$endo2), rownames(tmp))
    oupGEX$endo2 = cbind(oupGEX$endo2[commonGenes, ], tmp[commonGenes, ])
  }
}
gc()


## Kidney_Kuppe_4059315 (zenodo)
# Discard cells from CDpX libraries
tmpMeta = fread("rawData/Kidney_Kuppe_4059315/CD10negative/kidneyMap_UMI_counts_colData.txt")
tmpMeta$cell = paste0("kidney1-norm-", tmpMeta$`Patient ID`, "_cell", 1:nrow(tmpMeta))
tmpMeta$cell = gsub("-norm-CDm8",  "-CKD-CDm8", tmpMeta$cell)
tmpMeta$cell = gsub("-norm-CDm9",  "-CKD-CDm9", tmpMeta$cell)
tmpMeta$cell = gsub("-norm-CDm10", "-CKD-CDm10", tmpMeta$cell)
tmpMeta$cell = gsub("-norm-CDm11", "-CKD-CDm11", tmpMeta$cell)
tmpMeta$cell = gsub("-norm-CDm12", "-CKD-CDm12", tmpMeta$cell)
tmpMeta$cell = gsub("-norm-CDm13", "-CKD-CDm13", tmpMeta$cell)
tmpGene = fread("rawData/Kidney_Kuppe_4059315/CD10negative/kidneyMap_UMI_counts_rowData.txt")
tmp = readMM(file = "rawData/Kidney_Kuppe_4059315/CD10negative/kidneyMap_UMI_counts.mtx")
colnames(tmp) = tmpMeta$cell
rownames(tmp) = tmpGene$Gene.Symbol
tmpMeta = tmpMeta[grep("CDm", `Patient ID`)]
tmp = tmp[, tmpMeta$cell]
# unique(tstrsplit(tmpMeta$cell, "_")[[1]])
oupGEX$kidney1 = tmp
oupGEX$kidney1 = as(oupGEX$kidney1, "dgCMatrix")
for(i in unique(tstrsplit(colnames(oupGEX$kidney1), "_")[[1]])){
  print(paste0(i, " | ", dim(oupGEX$kidney1[, grep(i, colnames(oupGEX$kidney1))])[1], 
               " | ", dim(oupGEX$kidney1[, grep(i, colnames(oupGEX$kidney1))])[2]))
}


## Kidney_Lake_atlas.kpmp.org
inpFiles = LoadH5Seurat(
  "rawData/Kidney_Lake_atlas.kpmp.org/521c5b34-3dd0-4871-8064-61d3e3f1775a_PREMIERE_Alldatasets_08132021.h5Seurat")
inpFiles$newlabel = inpFiles$sampletype
inpFiles@meta.data[inpFiles$sampletype == "DKD", ]$newlabel = "CKD"
inpFiles@meta.data[inpFiles$sampletype == "HCKD", ]$newlabel = "CKD"
inpFiles@meta.data[inpFiles$sampletype == "LD", ]$newlabel = "norm"
inpFiles$newlabel = paste0("kidney2-", inpFiles$newlabel, "-", rownames(inpFiles@meta.data))
oupGEX$kidney2 = inpFiles@assays$RNA@counts
colnames(oupGEX$kidney2) = inpFiles$newlabel
oupGEX$kidney2 = drop0(oupGEX$kidney2)
for(i in unique(tstrsplit(colnames(oupGEX$kidney2), "_")[[1]])){
  print(paste0(i, " | ", dim(oupGEX$kidney2[, grep(i, colnames(oupGEX$kidney2))])[1], 
               " | ", dim(oupGEX$kidney2[, grep(i, colnames(oupGEX$kidney2))])[2]))
}
tmp = unique(data.table(sampleID = tstrsplit(inpFiles$newlabel, "_")[[1]],
                        status = tstrsplit(inpFiles$newlabel, "-")[[2]],
                        gender = inpFiles$gender,
                        age = inpFiles$age))



## Kidney_Malone_GSE145927
inpFiles = list.files("rawData/Kidney_Malone_GSE145927/", pattern = ".csv.gz")
oupGEX$kidney3 = NULL
for(i in inpFiles){
  tmp = fread(paste0("rawData/Kidney_Malone_GSE145927/", i))
  tmp = as(as.matrix(data.frame(tmp[, -1], row.names = tmp$V1)), "dgCMatrix")
  tmpID = gsub("rawcounts.csv.gz", "", strsplit(i, "_")[[1]][2])
  if(tmpID %in% c("day5", "day28")){
    tmpID = paste0("kidney3-norm-", tmpID)
  } else {
    tmpID = paste0("kidney3-rej-", tmpID)
  }
  print(paste0(tmpID, " | ", dim(tmp)[1], " | ", dim(tmp)[2]))
  colnames(tmp) = paste0(tmpID, "_", colnames(tmp))
  if(is.null(oupGEX$kidney3)){
    oupGEX$kidney3 = tmp
  } else {
    commonGenes = intersect(rownames(oupGEX$kidney3), rownames(tmp))
    oupGEX$kidney3 = cbind(oupGEX$kidney3[commonGenes, ], tmp[commonGenes, ])
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
  if(j == "skin1"){
    oupMeta = oupMeta[pctMT < 40 & pctHB < 0.1]
  } else {
    oupMeta = oupMeta[pctMT < 20 & pctHB < 0.1]
  }
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
  if(j == "kidney1"){
    tmp = data.table(p7$data)[features.plot == "CD68" & pct.exp > 10 &
                                avg.exp.scaled > 0]$id; print(tmp)
    seu = subset(seu, idents = tmp)
  } else {
    tmp = data.table(p7$data)[features.plot == "CD68" & pct.exp > 40 &
                                avg.exp.scaled > 0]$id; print(tmp)
    seu = subset(seu, idents = tmp)
  }
  oupGEX[[j]] = oupGEX[[j]][, colnames(seu)]
  oupGEX[[j]] = oupGEX[[j]][rowSums(oupGEX[[j]]) > 0, ]
  oupNum[study == j]$r2nCell = ncol(oupGEX[[j]])
  oupNum[study == j]$r2nGene = nrow(oupGEX[[j]])
  fwrite(oupNum, sep = "\t", file = "res_0prepData1/oupNUMr2b.txt")
  saveRDS(oupGEX[[j]], file = paste0("res_0prepData1rds/oupGEXr2_", j, ".rds"))
  gc()
}




