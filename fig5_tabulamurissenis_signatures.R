library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(ggrepel)
library(dplyr) 
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(rio)
library(UCell)
library(TabulaMurisData)
library(ExperimentHub)
library(ggforce)
library(rstatix)
library(SingleCellExperiment)


setwd("~/Documents/nBox/eLife_MAM/Aging Analysis/Mouse/")
mam.sig <-
  import("~/Documents/nBox/eLife_MAM/Objects/mam_DEsumm.txt")
rownames(mam.sig) <- mam.sig$gene
mam.pos <-
  unique(babelgene::orthologs(mam.sig[mam.sig$mamGenes, "gene"], species = "mouse", human = T)$symbol)
mam.neg <-
  unique(babelgene::orthologs(mam.sig[mam.sig$avgLFC < -0.5 &
                                        mam.sig$zDiffProp < -1.96, "gene"], species = "mouse", human = T)$symbol)
mam.metab <-
  unique(babelgene::orthologs(mam.sig[mam.sig$mamGenesMetab, "gene"], species = "mouse", human = T)$symbol)
mam.ecm <-
  unique(babelgene::orthologs(mam.sig[mam.sig$mamGenesECM, "gene"], species = "mouse", human = T)$symbol)

# sceasy::convertFormat(obj = "tabula-muris-senis-droplet-processed-official-annotations.h5ad", from = "anndata", to = "seurat", outFile = "tabula_muris_senis.RDS")
tms <- readRDS(file = "tabula_muris_senis.RDS")

# Lung ####
lung.tms <- subset(
  tms,
  tissue == "Lung" &
    cell_ontology_class %in% c(
      "classical monocyte",
      "intermediate monocyte",
      "alveolar macrophage",
      "lung macrophage",
      "non-classical monocyte"
    )
)
lung.tms <- NormalizeData(lung.tms)
lung.tms <- FindVariableFeatures(lung.tms)
lung.tms <- ScaleData(lung.tms)
lung.tms <- RunPCA(lung.tms)
lung.tms <- FindNeighbors(lung.tms, dims = 1:50)
lung.tms <- FindClusters(lung.tms, resolution = 0.5)
lung.tms <- RunUMAP(lung.tms, dims = 1:50)
DimPlot(lung.tms, label = T, raster = F) + NoLegend()

# lung.tms.markers <- presto::wilcoxauc(lung.tms)
lung.tms.markers <- FindAllMarkers(lung.tms, only.pos = T)
## All macrophages ####
lung.tms <-
  AddModuleScore(object = lung.tms,
                 features = list(mam.pos),
                 name = "MAMPos")
lung.tms <-
  AddModuleScore(object = lung.tms,
                 features = list(mam.metab),
                 name = "Metabolic")
lung.tms <-
  AddModuleScore(object = lung.tms,
                 features = list(mam.ecm),
                 name = "ECM")
lung.tms <-
  AddModuleScore(object = lung.tms,
                 features = list(mam.neg),
                 name = "MAMNeg")


lung.tms$age_group <-
  ifelse(lung.tms$age %in% c("1m", "3m"),
         "Young (1m, 3m)",
         "Aged (21m, 30m)")
mon.clusters.lung <-
  unique(lung.tms.markers[lung.tms.markers$gene %in% c("Fcn1", "S100a8", "S100a9", "Lst1", "Fcgr3a") &
                            lung.tms.markers$p_val_adj < 0.05 &
                            lung.tms.markers$avg_log2FC > 0.5, "cluster"])

scored.macrophages.lung <- subset(lung.tms, idents = c(mon.clusters.lung), invert = T)

vln.dfs.lung <- list()
for (signatures in c("MAMPos1",
                     "MAMNeg1",
                     "ECM1",
                     "Metabolic1")) {
  vln.dfs.lung[[signatures]] <-
    data.frame(
      cell = colnames(scored.macrophages.lung),
      signature = signatures,
      score = scored.macrophages.lung@meta.data[, signatures],
      age = scored.macrophages.lung$age_group
    )
}
vln.df.lung <- do.call(rbind, vln.dfs.lung)
vln.df.lung$age <-
  factor(vln.df.lung$age, levels = c("Young (1m, 3m)", "Aged (21m, 30m)"))


vln.df.lung$signature <-
  gsub("ECM1", "ECM", gsub("Metabolic1", "Metabolic", gsub(
    "MAMPos1", "SPP1+MAM+", gsub(
      "MAMNeg1",
      "SPP1+MAM-",
      vln.df.lung$signature)
    )
  ))
vln.df.lung$signature <- factor(vln.df.lung$signature, levels = c("SPP1+MAM-", "SPP1+MAM+", "ECM", "Metabolic"))
lung.violin <-
  ggplot(vln.df.lung, aes(x = age, y = score,  fill = age)) +
  geom_violin(
    adjust = 2.5,
    width = 0.8,
    scale = 'width',
    aes(fill = age),
    position = position_dodge(width = 0.9)
  ) +
  geom_point(
    position = position_jitterdodge(seed = 1, dodge.width = 0.9),
    color = "black",
    size = 0.5
  ) +
  ylim(0, 1) +
  scale_fill_manual(values = c(
    "Young (1m, 3m)" = "orange2",
    "Aged (21m, 30m)" = "purple1"
  )) +
  facet_wrap( ~ signature, nrow = 1, strip.position = "bottom") +
  stat_compare_means(
    comparisons = list(c("Aged (21m, 30m)", "Young (1m, 3m)")),
    method = "wilcox.test",
    method.args = list("alternative" = "two.sided"),
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.15,
    label.y = 0.8,
  ) +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 20,
    colour = "white",
    shape = 95
  ) +
  theme_classic(base_size = 30) +
  xlab("") +
  ylab("Signature score in Lung macrophages") +
  labs(fill = "Age") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 32),
    panel.border = element_rect(
      color = "white",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "white", size = 1)
  )

pdf(file = "plots/violin_plots.pdf", width = 25, height = 10)
lung.violin
hcdev.off()



# Kidney ####
kidney.tms <-
  subset(tms,
         tissue == "Kidney" & cell_ontology_class %in% c("macrophage"))

kidney.tms <- NormalizeData(kidney.tms)
kidney.tms <- FindVariableFeatures(kidney.tms)
kidney.tms <- ScaleData(kidney.tms)
kidney.tms <- RunPCA(kidney.tms)
kidney.tms <- FindNeighbors(kidney.tms, dims = 1:50)
kidney.tms <- FindClusters(kidney.tms, resolution = 0.5)
kidney.tms <- RunUMAP(kidney.tms, dims = 1:50)
DimPlot(kidney.tms, label = T, raster = F) + NoLegend()

# kidney.tms.markers <- presto::wilcoxauc(kidney.tms)
kidney.tms.markers <- FindAllMarkers(kidney.tms, only.pos = T)

## All Macrophages ####
kidney.tms <-
  AddModuleScore(object = kidney.tms,
                 features = list(mam.pos),
                 name = "MAMPos")
kidney.tms <-
  AddModuleScore(object = kidney.tms,
                 features = list(mam.metab),
                 name = "Metabolic")
kidney.tms <-
  AddModuleScore(object = kidney.tms,
                 features = list(mam.ecm),
                 name = "ECM")
kidney.tms <-
  AddModuleScore(object = kidney.tms,
                 features = list(mam.neg),
                 name = "MAMNeg")


kidney.tms$age_group <-
  ifelse(kidney.tms$age %in% c("1m", "3m"),
         "Young (1m, 3m)",
         "Aged (21m, 30m)")
mon.clusters.kidney <-
  unique(kidney.tms.markers[kidney.tms.markers$gene %in% c("Fcn1", "S100a8", "S100a9", "Lst1", "Fcgr3a") &
                              kidney.tms.markers$p_val_adj < 0.05 &
                              kidney.tms.markers$avg_log2FC > 0.5, "cluster"])

scored.macrophages.kidney <-
  subset(kidney.tms,
         idents = c(mon.clusters.kidney),
         invert = T)

vln.dfs.kidney <- list()
for (signatures in c("MAMPos1",
                     "MAMNeg1",
                     "ECM1",
                     "Metabolic1"
                     )) {
  vln.dfs.kidney[[signatures]] <-
    data.frame(
      cell = colnames(scored.macrophages.kidney),
      signature = signatures,
      score = scored.macrophages.kidney@meta.data[, signatures],
      age = scored.macrophages.kidney$age_group
    )
}
vln.df.kidney <- do.call(rbind, vln.dfs.kidney)
vln.df.kidney$age <-
  factor(vln.df.kidney$age, levels = c("Young (1m, 3m)", "Aged (21m, 30m)"))


vln.df.kidney$signature <-
  gsub("ECM1", "ECM", gsub("Metabolic1", "Metabolic", gsub(
    "MAMPos1", "SPP1+MAM+", gsub(
      "MAMNeg1",
      "SPP1+MAM-",
      vln.df.kidney$signature)
    )
  ))

vln.df.kidney$signature <- factor(vln.df.kidney$signature, levels = c("SPP1+MAM-", "SPP1+MAM+", "ECM", "Metabolic"))
kidney.violin <-
  ggplot(vln.df.kidney, aes(x = age, y = score,  fill = age)) +
  geom_violin(
    adjust = 2.5,
    width = 0.8,
    scale = 'width',
    aes(fill = age),
    position = position_dodge(width = 0.9)
  ) +
  geom_point(
    position = position_jitterdodge(seed = 1, dodge.width = 0.9),
    color = "black",
    size = 0.5
  ) +
  ylim(0, 1) +
  scale_fill_manual(values = c(
    "Young (1m, 3m)" = "orange2",
    "Aged (21m, 30m)" = "purple1"
  )) +
  facet_wrap( ~ signature, nrow = 1, strip.position = "bottom") +
  stat_compare_means(
    comparisons = list(c("Aged (21m, 30m)", "Young (1m, 3m)")),
    method = "wilcox.test",
    method.args = list("alternative" = "two.sided"),
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.15,
    label.y = 0.8,
  ) +
  stat_summary(
    fun.y = median,
    geom = "point",
    size = 20,
    colour = "white",
    shape = 95
  ) +
  theme_classic(base_size = 30) +
  xlab("") +
  ylab("Signature score in kidney macrophages") +
  labs(fill = "Age") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 32),
    panel.border = element_rect(
      color = "white",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "white", size = 1)
  )  

pdf(file = "plots/violin_plots_kidney.pdf", width = 25, height = 10)
kidney.violin
dev.off()



