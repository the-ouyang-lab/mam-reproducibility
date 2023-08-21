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

setwd("~/Documents/nBox/eLife_MAM/Aging Analysis/")

# Signatures #####
mam.sig <- import("../Objects/mam_DEsumm.txt")
rownames(mam.sig) <- mam.sig$gene
mam.pos <- mam.sig[mam.sig$mamGenes, "gene"]
mam.neg <-
  mam.sig[mam.sig$avgLFC < -0.5 & mam.sig$zDiffProp < -1.96, "gene"]
mam.metab <- mam.sig[mam.sig$mamGenesMetab, "gene"]
mam.ecm <- mam.sig[mam.sig$mamGenesECM, "gene"]
homeo <- import("../Objects/homeo_DEsumm-updated.txt")
rownames(homeo) <- homeo$gene
homeo.sig <- homeo[homeo$homeoGenes == T, "gene"]

# Reading in dataset ####
hgnc <- import("../../Kidney Project/hgnc_genes_list_180122.txt")
hgnc <-
  hgnc[hgnc$Status == "Approved"  &
         length(hgnc$`Ensembl ID(supplied by Ensembl)`) > 0,]

hlca <- readRDS("Human Lung Cell Atlas/hlca_core.rds")


# Filtering and preprocessing macrophages #####
macs.mons.hlca <-
  subset(
    hlca,
    cell_type %in% c(
      "alveolar macrophage",
      "classical monocyte",
      "elicited macrophage",
      "non-classical monocyte",
      "lung macrophage"
    )
  )

rownames(macs.mons.hlca@assays$RNA@counts) <- plyr::mapvalues(
  x = rownames(macs.mons.hlca@assays$RNA@counts),
  from = hgnc$`Ensembl ID(supplied by Ensembl)`,
  to = hgnc$`Approved symbol`,
  warn_missing = F
)
rownames(macs.mons.hlca@assays$RNA@data) <- plyr::mapvalues(
  x = rownames(macs.mons.hlca@assays$RNA@counts),
  from = hgnc$`Ensembl ID(supplied by Ensembl)`,
  to = hgnc$`Approved symbol`,
  warn_missing = F
)


macs.mons.hlca <- NormalizeData(macs.mons.hlca)
macs.mons.hlca <- FindVariableFeatures(macs.mons.hlca)
macs.mons.hlca <- ScaleData(macs.mons.hlca)
macs.mons.hlca <- RunPCA(macs.mons.hlca)
macs.mons.hlca <- FindNeighbors(macs.mons.hlca, dims = 1:50)
macs.mons.hlca <- FindClusters(macs.mons.hlca, resolution = 0.5)
macs.mons.hlca <- RunUMAP(macs.mons.hlca, dims = 1:50)
DimPlot(macs.mons.hlca, label = T, raster = F) + NoLegend()

macs.mons.hlca.markers <- presto::wilcoxauc(macs.mons.hlca, only.pos = T)
# macs.mons.hlca.markers <- FindAllMarkers(macs.mons.hlca, only.pos = T)

monocyte.clusters <-
  unique(macs.mons.hlca.markers[macs.mons.hlca.markers$feature %in% 
                             c("FCN1", "S100A8", "S100A9", "LST1", "FCGR3A") &
                             macs.mons.hlca.markers$padj < 0.05 &
                             macs.mons.hlca.markers$logFC > 0.5, "group"])

# Scoring macrophages for the signatures ####
macs.mons.hlca <-
  AddModuleScore(object = macs.mons.hlca,
                 features = list(mam.pos),
                 name = "MAMPos")
macs.mons.hlca <-
  AddModuleScore(object = macs.mons.hlca,
                 features = list(mam.metab),
                 name = "MAMMet")
macs.mons.hlca <-
  AddModuleScore(object = macs.mons.hlca,
                 features = list(mam.ecm),
                 name = "MAMECM")
macs.mons.hlca <-
  AddModuleScore(object = macs.mons.hlca,
                 features = list(mam.neg),
                 name = "MAMNeg")
macs.mons.hlca <-
  AddModuleScore(object = macs.mons.hlca,
                 features = list(homeo.sig),
                 name = "Homeostatic")

donors.to.keep <-
  names(table(macs.mons.hlca$donor_id))[table(macs.mons.hlca$donor_id) > 10]

scored.lung.macrophages <-
  subset(macs.mons.hlca,
         donor_id %in% donors.to.keep)

scored.lung.macrophages <-
  subset(scored.lung.macrophages,
         idents = c(monocyte.clusters),
         invert = T)

pseudo.lung <- data.frame(
  row.names = colnames(scored.lung.macrophages),
  MAM.pos = scored.lung.macrophages$MAMPos1,
  MAM.met = scored.lung.macrophages$MAMMet1,
  MAM.ecm = scored.lung.macrophages$MAMECM1,
  MAM.neg = scored.lung.macrophages$MAMNeg1,
  homeostatic = scored.lung.macrophages$Homeostatic1,
  Age = as.numeric(scored.lung.macrophages$age_or_mean_of_age_range),
  Smoking = scored.lung.macrophages$smoking_status,
  Sex = scored.lung.macrophages$sex,
  Patient = scored.lung.macrophages$donor_id
)

pseudo.lung$Bins <- cut(pseudo.lung$Age, c(20, 40, 60, 80))
pseudo.lung$MAM.pos <- scale(pseudo.lung$MAM.pos)
pseudo.lung$MAM.neg <- scale(pseudo.lung$MAM.neg)
pseudo.lung$homeostatic <- scale(pseudo.lung$homeostatic)
pseudo.lung$MAM.ecm <- scale(pseudo.lung$MAM.ecm)
pseudo.lung$MAM.met <- scale(pseudo.lung$MAM.met)


pseudo.median.lung <-
  na.omit(unique(
    pseudo.lung %>%
      group_by(Patient) %>%
      reframe(
        median.MAM.pos = median(MAM.pos, na.rm = TRUE),
        median.MAM.met = median(MAM.met, na.rm = TRUE),
        median.MAM.ecm = median(MAM.ecm, na.rm = TRUE),
        median.MAM.neg = median(MAM.neg, na.rm = TRUE),
        median.homeo = median(homeostatic, na.rm = TRUE),
        Sex,
        Smoking,
        Age,
        Bins
      )
  ))
pseudo.median.lung$Smoking <- stringr::str_to_sentence(pseudo.median.lung$Smoking)
# Lineplots ####
p1 <- ggscatter(
  data = pseudo.median.lung,
  x = "Age",
  y = "median.homeo",
  cor.coef = TRUE,
  add = "reg.line",
  conf.int = TRUE,
  # facet.by = "Smoking",
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coeff.args = list(
    method = "spearman",
    size = 10,
    label.y = 3,
    label.sep = ", "
  )
) +
  xlim(20, 80) +
  ylim(-3, 3) +
  theme_classic(base_size = 30) +
  theme(axis.text = element_text(size = 35)) +
  xlab("") +
  ylab("Normalized Pseudobulk score") +
  ggtitle("Homeostatic RNASE1+")
p1
p2 <- ggscatter(
  data = pseudo.median.lung,
  x = "Age",
  y = "median.MAM.neg",
  cor.coef = TRUE,
  add = "reg.line",
  # facet.by = "Smoking",
  conf.int = TRUE,
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coeff.args = list(
    method = "spearman",
    size = 10,
    label.y = 3,
    label.sep = ", "
  )
) +
  xlim(20, 80) +
  ylim(-3, 3) +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  xlab("") +
  ylab("") +
  ggtitle("SPP1+MAM-")

p3 <- ggscatter(
  data = pseudo.median.lung,
  x = "Age",
  y = "median.MAM.pos",
  cor.coef = TRUE,
  add = "reg.line",
  conf.int = TRUE,
  # facet.by = "Smoking",
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coeff.args = list(
    method = "spearman",
    size = 10,
    label.y = 3,
    label.sep = ", "
  )
)  +
  xlim(20, 80) +
  ylim(-3, 3) +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  xlab("Age") +
  ylab("") +
  ggtitle("SPP1+MAM+")

p4 <- ggscatter(
  data = pseudo.median.lung,
  x = "Age",
  y = "median.MAM.ecm",
  cor.coef = TRUE,
  add = "reg.line",
  conf.int = TRUE,
  # facet.by = "Smoking",
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coeff.args = list(
    method = "spearman",
    size = 10,
    label.y = 3,
    label.sep = ", "
  )
)  +
  xlim(20, 80) +
  ylim(-3, 3) +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  xlab("") +
  ylab("") +
  ggtitle("ECM")

p5 <- ggscatter(
  data = pseudo.median.lung,
  x = "Age",
  y = "median.MAM.met",
  cor.coef = TRUE,
  add = "reg.line",
  conf.int = TRUE,
  # facet.by = "Smoking",
  add.params = list(color = "blue", fill = "lightgray"),
  cor.coeff.args = list(
    method = "spearman",
    size = 10,
    label.y = 3,
    label.sep = ", "
  )
)  +
  xlim(20, 80) +
  ylim(-3, 3) +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  xlab("") +
  ylab("") +
  ggtitle("Metabolic")

pdf(file = "Human Lung Cell Atlas/scatter_plots.pdf",
    width = 35,
    height = 8)
ggarrange(plotlist = list(p1, p2, p3, p4, p5), nrow = 1)
dev.off()


# Boxplots ####
p6 <- ggplot(pseudo.median.lung, aes(x = Bins, y = median.homeo)) +
  geom_boxplot(aes(fill = Bins),) +
  geom_point(position = position_jitter(width = 0.2)) +
  ggtitle("Homeostatic RNASE1+") +
  theme_classic(base_size = 30) +
  theme(axis.text = element_text(size = 35)) +
  scale_x_discrete(labels = c("20-40", "40-60", "60-80")) +
  xlab("") + ylab("Normalized Pseudobulk score")  +
  ylim(-5, 6) +
  stat_compare_means(
    comparisons = list(c("(40,60]", "(20,40]"), c("(60,80]", "(20,40]")),
    method = "wilcox.test",
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.12,
    label.y = 4
  ) + NoLegend()

p7 <-
  ggplot(pseudo.median.lung, aes(x = Bins, y = median.MAM.neg)) +
  geom_boxplot(aes(fill = Bins)) +
  geom_point(position = position_jitter(width = 0.2)) +
  ggtitle("SPP1+MAM-") +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  scale_x_discrete(labels = c("20-40", "40-60", "60-80")) +
  xlab("") + ylab("")  +
  ylim(-5, 6) +
  stat_compare_means(
    comparisons = list(c("(40,60]", "(20,40]"), c("(60,80]", "(20,40]")),
    method = "wilcox.test",
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.15,
    label.y = 4
  ) + NoLegend()

p8 <-
  ggplot(pseudo.median.lung, aes(x = Bins, y = median.MAM.pos)) +
  geom_boxplot(aes(fill = Bins)) +
  geom_point(position = position_jitter(width = 0.2)) +
  ggtitle("SPP1+MAM+") +
  theme_classic(base_size = 30) +
  theme(
    axis.text.y = element_blank(),
    axis.text = element_text(size = 35),
    axis.title.x = element_text(size = 35)
  ) +
  scale_x_discrete(labels = c("20-40", "40-60", "60-80")) +
  xlab("Age") + ylab("")  +
  ylim(-5, 6) +
  stat_compare_means(
    comparisons = list(c("(40,60]", "(20,40]"), c("(60,80]", "(20,40]")),
    method = "wilcox.test",
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.20,
    label.y = 4
  ) + NoLegend()

p9 <-
  ggplot(pseudo.median.lung, aes(x = Bins, y = median.MAM.ecm)) +
  geom_boxplot(aes(fill = Bins)) +
  geom_point(position = position_jitter(width = 0.2)) +
  ggtitle("ECM") +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  scale_x_discrete(labels = c("20-40", "40-60", "60-80")) +
  xlab("") + ylab("")  +
  ylim(-5, 6) +
  stat_compare_means(
    comparisons = list(c("(40,60]", "(20,40]"), c("(60,80]", "(20,40]")),
    method = "wilcox.test",
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.25,
    label.y = 4
  ) + NoLegend()

p10 <-
  ggplot(pseudo.median.lung, aes(x = Bins, y = median.MAM.met)) +
  geom_boxplot(aes(fill = Bins)) +
  geom_point(position = position_jitter(width = 0.2)) +
  ggtitle("Metabolic") +
  theme_classic(base_size = 30) +
  theme(axis.text.y = element_blank(),
        axis.text = element_text(size = 35)) +
  scale_x_discrete(labels = c("20-40", "40-60", "60-80")) +
  xlab("") + ylab("")  +
  ylim(-5, 6) +
  stat_compare_means(
    comparisons = list(c("(40,60]", "(20,40]"), c("(60,80]", "(20,40]")),
    method = "wilcox.test",
    size = 10,
    bracket.size = 0.1,
    tip.length = 0,
    step.increase = 0.25,
    label.y = 4
  ) + NoLegend()

pdf(file = "Human Lung Cell Atlas/box_plots.pdf",
    width = 36,
    height = 10)
ggarrange(plotlist = list(p6, p7, p8, p9, p10), nrow = 1)
dev.off()


# Regression ####

lm.homeo <-
  lm(median.homeo ~ Age + Sex + Smoking, data = pseudo.median.lung)
summary(lm.homeo)

lm.mam.neg <-
  lm(median.MAM.neg ~ Age + Sex + Smoking, data = pseudo.median.lung)
summary(lm.mam.neg)

lm.mam.pos <-
  lm(median.MAM.pos ~ Age + Sex + Smoking, data = pseudo.median.lung)
summary(lm.mam.pos)

lm.mam.ecm <-
  lm(median.MAM.ecm ~ Age + Sex + Smoking, data = pseudo.median.lung)
summary(lm.mam.ecm)

lm.mam.met <-
  lm(median.MAM.met ~ Age + Sex + Smoking, data = pseudo.median.lung)
summary(lm.mam.met)
