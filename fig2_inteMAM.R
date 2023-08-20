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
if(!dir.exists("res_1inteData/mergeClust/")){dir.create("res_1inteData/mergeClust/")}

hmapCol <- colorRampPalette(c("blue", "white", "red"))(21)
hmapBrk <- seq(-1, 1, length.out = 21)
genes.to.label <- c(
  'SELENOP','C1QA','RNASE1' ,'C1QC','FCN1', 'SDC2', 'PLA2G7','S100A9','S100A8',
  'VCAN', 'APOE','GPNMB','CD9','STAB1','SPP1', 'FABP5', 'TREM2', 'F13A1', 
  'CCL20', 'FOLR2', 'INHBA',  'CHI3L1',  'MATK','PTGS2','IL1B', 'FABP4', 'LYVE1', 
  'APOC1', 'CCL4','S100A12', 'CTSK', 'LST1', 'CDKN1C', 'FCGR3A',
  'CCL3', 'DAB2','CCL4L2', 'PHLDA3', 'EGR1', 'ATF3', 'VEGFA', 'SEMA6B', 
  'IFI44L', 'IFIT3', 'IFIT1', 'MT1G', 'MT1X', 'MT1H')

# liver   10401
# lung   184297
# heart   24461
# skin     2776
# endo     6030
# kidney   7965



### 1. Start Integration
oupSPP1gex = list()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  # Read in seurat object
  seu = readRDS(paste0("res_1inteData/", i, ".rds"))
  # Subset SPP1 cells for next step
  tmp = levels(Idents(seu))
  print(tmp[grep("spp1-", tmp)])
  seu = subset(seu, idents = tmp[grep("spp1-", tmp)])
  oupSPP1gex[[i]] = seu@assays$RNA@counts
  oupSPP1gex[[i]] = oupSPP1gex[[i]][rowSums(oupSPP1gex[[i]]) > 0, ]
}

### Integrate SPP1 potato
# Get genes that are expressed in 2 or more datasets
commonGenes = c()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  commonGenes = c(commonGenes, rownames(oupSPP1gex[[i]]))}
commonGenes = names(table(commonGenes))[table(commonGenes) >= 6]
for(i in c("liver","lung","heart","skin","endo","kidney")){
  tmpNewG = setdiff(commonGenes, rownames(oupSPP1gex[[i]]))
  if(length(tmpNewG) > 0){
    tmpNewM = matrix(0, ncol = ncol(oupSPP1gex[[i]]), nrow = length(tmpNewG))
    rownames(tmpNewM) = tmpNewG; colnames(tmpNewM) = colnames(oupSPP1gex[[i]])
    oupSPP1gex[[i]] = rbind(oupSPP1gex[[i]], tmpNewM)
  }
  oupSPP1gex[[i]] = oupSPP1gex[[i]][commonGenes, ]
}

# Start SCTransform
oupSeuList = list()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  tmp = oupSPP1gex[[i]]
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
seu <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = (4:20)/20, verbose = FALSE)

# Plot 
colIdent = brewer.pal(6, "Dark2")
ggOut = list()
ggOut[["sample"]] = DimPlot(seu, group.by = "tissue", cols = colIdent, 
                            shuffle = T, pt.size = 0.25) + NoLegend() + coord_fixed() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())
# tmpG = c("SPP1", "MMP7", "MMP9", "TREM2", "FCN1")
tmpG = c("SPP1", "TREM2", "MMP9", 
         "GSN", "MATK", "MFSD12", "PKM", "ENO1")
for(j in tmpG){
  if(j %in% rownames(seu)){DefaultAssay(seu) = "integrated"}
  if(!j %in% rownames(seu)){DefaultAssay(seu) = "RNA"}
  DefaultAssay(seu) = "SCT"
  ggOut[[j]] = FeaturePlot(seu, features = j, order = F, pt.size = 0.1) + 
    scale_color_gradientn(colors = c(
      "grey85", brewer.pal(9, "Reds")[c(1:9,9,9)])) + coord_fixed() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank())
}
ggsave(wrap_plots(ggOut), height = 12, width = 16,
       filename = paste0("res_1inteData/SPP1_umap.png"))
saveRDS(seu, file = paste0("res_1inteData/SPP1seu.rds"))



### 2. Identify MAMs
## A. Calculate silouette score and determine best res
oupSil = data.table()
tmpD = dist(x = Embeddings(seu)[, 1:30])
for(j in (4:20)/20){
  tmp = silhouette(x = as.numeric(seu@meta.data[[
    paste0("integrated_snn_res.", j)]]), dist = tmpD)
  tmp = data.table(tissue = "SPP1potato", resolution = j,
                   silScore = mean(summary(tmp)$clus.avg.widths))
  oupSil = rbindlist(list(oupSil, tmp))
}
tmpD = oupSil[tissue == "SPP1potato" & resolution >= 0.3][which.max(silScore)]
p1 = ggplot(oupSil[tissue == "SPP1potato"],
            aes(resolution, silScore, label = round(silScore, 3))) +
  geom_point() + geom_text(vjust = 1.35, size = 4) +
  geom_point(data = tmpD, shape = 1, color = "red", size = 6) +
  ylab("silhouette score") + theme_classic(base_size = 20)
ggsave(p1, height = 6, width = 6,
       filename = paste0("res_1inteData/mergeClust/SPP1_silScore.png"))

## B. Merge clusters
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
pheatmap::pheatmap(tmpCor, breaks = hmapBrk, color = hmapCol, border_color = "white", 
                   display_numbers = T, number_format = "%.02f", number_color = "black",
                   fontsize_row = 14, fontsize_col = 14,
                   fontsize_number = 14, height = 8, width = 8, angle_col = 45, 
                   filename = paste0("res_1inteData/mergeClust/SPP1_spearman.png"))

## C. Do some plots
seu$clust = paste0("C", as.character(seu$seurat_clusters))
seu@meta.data[seu$clust == "C3", "clust"] = "C3-MAM"
seu$clust = factor(seu$clust, levels = c(
  "C3-MAM", paste0("C", 0:2), paste0("C", 4:8)))
seu$mam = "MAM"
seu@meta.data[seu$clust != "C3-MAM", "mam"] = "nonMAM"
seu$mam = factor(seu$mam, levels = c("MAM", "nonMAM"))
seu$tissue_mam = paste0(seu$tissue, "_", seu$mam)

Idents(seu) <- seu$clust
p1 = DimPlot(seu, group.by = "clust", label = TRUE, raster = FALSE,
             cols = colorRampPalette(brewer.pal(10, "Paired"))(uniqueN(seu$clust)))
p2 = DotPlot(seu, group.by = "clust", assay = "SCT",
             features = c("SPP1", "TREM2", "MMP9", 
                          "GSN", "MATK", "MFSD12", "PKM", "ENO1")) + 
  scale_color_gradientn(colors = c("grey85", brewer.pal(9, "Reds")))
ggsave(p1 + p2, height = 6, width = 14, 
       filename = paste0("res_1inteData/mergeClust/SPP1_clustGEX.png"))

ggData <- data.frame(prop.table(table(seu$tissue, seu$clust), margin = 2))
colnames(ggData) = c("tissue", "cluster", "value")
ggData$tissue = factor(ggData$tissue, levels = c(
  "liver","lung","heart","skin","endo","kidney"))
p1 <- ggplot(ggData, aes(cluster, value, fill = tissue)) +
  geom_col() + xlab("Cluster") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = brewer.pal(6, "Dark2")) + coord_flip() +
  theme_classic(base_size = 16)
ggsave(p1, height = 8, width = 6,
       filename = paste0("res_1inteData/mergeClust/SPP1_prop.png"))

## D. Derive MAM signature
oupDE = data.table()
# Subset then do DE
for(i in c("liver","lung","heart","skin","endo","kidney")){
  seuDE <- subset(seu, tissue == i)
  DefaultAssay(seuDE) <- "SCT"
  tmpMarkers <- FoldChange(seuDE, ident.1 = "C3-MAM", pseudocount.use = 0.1)
  tmpMarkers <- data.table(tissue = i, gene = rownames(tmpMarkers), tmpMarkers)
  tmpMarkers$diffProp = tmpMarkers$pct.1 - tmpMarkers$pct.2
  oupDE = rbindlist(list(oupDE, tmpMarkers))
}

# Calculate avgLFC, avgDiffProp, zDiffProp
oupDEsumm = oupDE[, .(avgLFC = sum(avg_log2FC) / 6,
                      avgDiffProp = sum(diffProp) / 6), by = "gene"]
oupDEsumm$zDiffProp = (oupDEsumm$avgDiffProp - mean(oupDEsumm$avgDiffProp)) /
  sd(oupDEsumm$avgDiffProp)
oupDEsumm$mamGenes = FALSE
oupDEsumm[zDiffProp > 1.96][avgLFC > 0.5]$mamGenes = TRUE
oupDEsumm = oupDEsumm[order(-mamGenes,-avgLFC)]
oupDEsumm[order(-avgLFC)][mamGenes == TRUE][1:20]

## E. Do enrichment of MAM signature
inpGS = data.table(msigdbr(species = "human", category = "C2"))
inpGS = inpGS[gs_subcat %in% c(
  "CP:BIOCARTA","CP:KEGG","CP:PID","CP:REACTOME"),
  c("gs_name","gene_symbol")]
inpLFC = oupDEsumm$avgLFC; names(inpLFC) = oupDEsumm$gene
inpLFC = sort(inpLFC, decreasing = T)
oupGSEA = clusterProfiler::GSEA(inpLFC, TERM2GENE = inpGS, seed = 42)
oupGSEA = data.table(oupGSEA@result)[order(-NES)]

# Plot top enrichments (NES > 1.9)
ggData = oupGSEA[NES > 1.8, -c(1,10,11)]
ggData$Description = gsub("REACTOME_|PID_|KEGG_", "", ggData$Description)
ggData$Description = gsub("_", " ", ggData$Description)
ggData$Description = stringr::str_to_title(ggData$Description)
ggData$mlog10Padj = -log10(ggData$p.adjust)
ggData$type = "nil"
ggData[grep("Extracellular|Collagen|Ecm|Integrin", Description)]$type = "ECM"
ggData[grep("Metabolism|Ppar|Insulin|Phosphorylation", Description)]$type = "Metab"
ggData[type == "nil"]$Description = ""
p1 <- ggplot(ggData, aes(NES, mlog10Padj, label = Description, 
                         size = setSize, color = type)) + 
  geom_point(alpha = 0.5) + geom_text_repel(size = 3, nudge_y = 0.05, force = 5) + 
  theme_classic(base_size = 18) + ylab("-log10(p adj.)") + xlim(c(1.75,2.15)) + 
  scale_size_continuous(range = c(5,15)) + ggtitle("Top 20 enriched terms") + 
  scale_color_manual(values = c("firebrick","dodgerblue","black"))
ggsave(p1, height = 6, width = 7,
       filename = paste0("res_1inteData/mam_gsea.png"))

# Extract ECM and metab genes and append to oupDEsumm
ggData = oupGSEA[NES > 1.6, -c(1,10)]
ggData$type = "nil"
ggData[grep("EXTRACELLULAR|COLLAGEN|ECM|INTEGRIN", Description)]$type = "ECM"
ggData[grep("METABOLISM|PPAR|UPAR|INSULIN|OXIDATIVE_PHOS|GLYCOLYSIS", Description)]$type = "Metab"
fwrite(oupGSEA,   sep = "\t", file = "res_1inteData/SPP1_FuncGSEA.txt")
fwrite(ggData,    sep = "\t", file = "res_1inteData/SPP1_FuncGSEAextra.txt")

tmpG = unique(unlist(strsplit(ggData[type == "ECM"]$core_enrichment, "\\/")))
oupDEsumm$mamGenesECM = FALSE
oupDEsumm[mamGenes == TRUE & gene %in% tmpG]$mamGenesECM = TRUE

tmpG = unique(unlist(strsplit(ggData[type == "Metab"]$core_enrichment, "\\/")))
oupDEsumm$mamGenesMetab = FALSE
oupDEsumm[mamGenes == TRUE & gene %in% tmpG]$mamGenesMetab = TRUE

## F. NABA database enrichment
inpGS = data.table(msigdbr(species = "human", category = "C2"))
inpGS = inpGS[grep("^NABA_", gs_name),
              c("gs_name","gene_symbol")]
inpGS = inpGS[!gs_name %in% c("NABA_CORE_MATRISOME",
                              "NABA_MATRISOME",
                              "NABA_MATRISOME_ASSOCIATED")]
oupNABA = clusterProfiler::enricher(
  oupDEsumm[mamGenes == TRUE]$gene, pvalueCutoff = 1,
  universe = oupDEsumm$gene, TERM2GENE = inpGS)
oupNABA = data.table(group = "mamGenes", oupNABA@result)

tmp = clusterProfiler::enricher(
  oupDEsumm[zDiffProp < -1.96][avgLFC < -0.5]$gene, pvalueCutoff = 1,
  universe = oupDEsumm$gene, TERM2GENE = inpGS)
tmp = data.table(group = "nonMam", tmp@result)
oupNABA = rbindlist(list(oupNABA, tmp))
oupNABA$mlog10P = -log10(oupNABA$pvalue)
oupNABA$Description = gsub("NABA_", "", oupNABA$Description)
oupNABA$Description = gsub("_", "\n", oupNABA$Description)
oupNABA$Description = stringr::str_to_title(oupNABA$Description)
oupNABA$Description = factor(oupNABA$Description, 
                             levels = rev(sort(unique(oupNABA$Description))))
fwrite(oupNABA,   sep = "\t", file = "res_1inteData/mam_FuncNABA.txt")

p1 <- ggplot(oupNABA[group == "mamGenes"], aes(mlog10P, Description)) + 
  geom_col() + ylab("") + xlab("-log10(p-value)") + xlim(c(0,1.8)) + 
  scale_y_discrete(limits = levels(oupNABA$Description)) + 
  theme_classic(base_size = 18) + ggtitle("Upregulated in SPP1+MAM+")
p2 <- ggplot(oupNABA[group == "nonMam"], aes(mlog10P, Description)) + 
  geom_col() + ylab("") + xlab("-log10(p-value)") + xlim(c(0,1.8)) + 
  scale_y_discrete(limits = levels(oupNABA$Description)) + 
  theme_classic(base_size = 18) + ggtitle("Upregulated in SPP1+MAM-")
ggsave(p1 + p2, height = 6, width = 12,
       filename = paste0("res_1inteData/mam_naba.png"))
# intersect(inpGS[gs_name == "NABA_SECRETED_FACTORS"]$gene_symbol, oupDEsumm[zDiffProp < -1.96][avgLFC < -0.5]$gene)
# intersect(inpGS[gs_name == "NABA_ECM_REGULATORS"]$gene_symbol, oupDEsumm[zDiffProp > 1.96][avgLFC > 0.5]$gene)
# ggData = oupGSEA[NES < -1.6, -c(1,10)]


### 2. Get other signatures
## Calculate marker gene
seu = PrepSCTFindMarkers(seu)
tmp = FindAllMarkers(seu, only.pos = T, logfc.threshold = 1)
tmp = data.table(tmp); tmp[, rank := frank(-avg_log2FC), by = "cluster"]
fwrite(tmp, sep = "\t", file = "res_1inteData/SPP1_markers.txt")


## Homeostatic signature
oupHomeo = data.table()
for(i in c("liver","lung","heart","skin","endo","kidney")){
  # Read in seurat object
  tmpSeu = readRDS(paste0("res_1inteData/", i, ".rds"))
  # Find markers for homeostatic
  homeoClust = levels(Idents(tmpSeu))[grep("^rnase", levels(Idents(tmpSeu)))]
  if(i == "endo"){homeoClust = "spp1-rnase-2_9"}
  DefaultAssay(tmpSeu) <- "SCT"
  tmpMarkers <- FoldChange(tmpSeu, ident.1 = homeoClust, pseudocount.use = 0.1)
  tmpMarkers <- data.table(tissue = i, gene = rownames(tmpMarkers), tmpMarkers)
  tmpMarkers$diffProp = tmpMarkers$pct.1 - tmpMarkers$pct.2
  oupHomeo = rbindlist(list(oupHomeo, tmpMarkers))
}

# Calculate avgLFC, avgDiffProp, zDiffProp
oupHOMsumm = oupHomeo[, .(avgLFC = sum(avg_log2FC) / 6,
                          avgDiffProp = sum(diffProp) / 6), by = "gene"]
oupHOMsumm$zDiffProp = (oupHOMsumm$avgDiffProp - mean(oupHOMsumm$avgDiffProp)) /
  sd(oupHOMsumm$avgDiffProp)
oupHOMsumm$homeoGenes = FALSE
oupHOMsumm[zDiffProp > 1.96][avgLFC > 0.5]$homeoGenes = TRUE
oupHOMsumm = oupHOMsumm[order(-homeoGenes,-avgLFC)]
oupHOMsumm[order(-avgLFC)][homeoGenes == TRUE][1:20]



### Final IO
fwrite(oupSil,    sep = "\t", file = "res_1inteData/mam_silScore.txt")
fwrite(oupDE,     sep = "\t", file = "res_1inteData/mam_DEres.txt")
fwrite(oupDEsumm, sep = "\t", file = "res_1inteData/mam_DEsumm.txt") # mam sig (158)
fwrite(oupHomeo,   sep = "\t", file = "res_1inteData/homeo_DEres.txt") # mam sig (99)
fwrite(oupHOMsumm, sep = "\t", file = "res_1inteData/homeo_DEsumm.txt") # mam sig (99)

saveRDS(oupSeuList, file = paste0("res_1inteData/SPP1seuList.rds"))
saveRDS(seu, file = paste0("res_1inteData/SPP1seu.rds"))




