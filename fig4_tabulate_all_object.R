setwd('~/project/Kevin/fig4/')

library(Seurat)
library(ggplot2)
library(patchwork)

SPP1seu <- readRDS('./data/SPP1seu.rds')

skin <- readRDS('./data/skin.rds')
kidney <- readRDS( './data/kidney.rds')
endo <- readRDS('./data/endo.rds')
lung <- readRDS('./data/lung.rds')
heart <- readRDS('./data/heart.rds')
liver <- readRDS('./data/liver.rds')

# Organizing meta-------------------
## Take for Kunal's code to re-annotate cells
#####SKIN############
skin <- RenameIdents(
  skin,
  "0" = "Trans",
  "3" = "Trans",
  "fcn1-2" = "FCN1+",
  "rnase-1" = "Homeostatic",
  "spp1-4" = "SPP1+MAM-"
)

skin <-
  SetIdent(skin, cells = Cells(subset(SPP1seu, tissue == "skin" &
                                        mam == "MAM")), value = "SPP1+MAM+")
skin$newclust <- Idents(skin)
skin$tissue <- 'skin'

table(skin$newclust)
# SPP1+MAM+       Trans       FCN1+ Homeostatic   SPP1+MAM-
#   8        1196         714         743         115

##### KIDNEY###################
kidney <- RenameIdents(
  kidney,
  "1_6" = "Trans",
  "2" = "Trans",
  "fcn1-0" = "FCN1+",
  "spp1-3" = "SPP1+MAM-",
  "spp1-5" = "SPP1+MAM-",
  "rnase-4" = "Homeostatic"
)

kidney <-
  SetIdent(kidney, cells = Cells(subset(SPP1seu, tissue == "kidney" &
                                          mam == "MAM")), value = "SPP1+MAM+")
kidney$newclust <- Idents(kidney)
kidney$tissue <- 'kidney'

table(kidney$newclust)

# SPP1+MAM+       Trans       FCN1+   SPP1+MAM- Homeostatic
# 103        3256        2138        1616         852

####ENDO##########################
# Get homeostatic cells
DefaultAssay(endo) <- 'SCT'
p1 <- FeaturePlot(endo, pt.size = 1, order = T,
                  features = c('FOLR2','LYVE1','MRC1', 'SPP1'))
p2 <- DimPlot(endo, group.by = 'clust', label = T, pt.size = 0.8) + NoLegend()
p1+p2+plot_layout(ncol = 3)

ggsave(p1+p2+plot_layout(ncol = 3), 
       filename='./output2/endo_homeo_marker.pdf',
       width = 14, height = 8)

# furthur split spp1-rnase-2_9
DefaultAssay(endo) <- 'integrated'
DimPlot(endo, group.by = 'integrated_snn_res.0.2', label = T)
Idents(endo) <- 'clust'
endo <- FindSubCluster(endo, cluster = "spp1-rnase-2_9", resolution = 0.1,
                       graph.name = 'integrated_snn', subcluster.name = "sub.cluster")
p3 <- DimPlot(endo, group.by = 'sub.cluster', label = T, pt.size = 0.8) + NoLegend()
p2+p3

ggsave(p2+p3+plot_layout(ncol = 2), 
       filename='./output2/endo_homeo_cluster.pdf',
       width = 10, height = 5)

## rename cell clusters
Idents(endo) <- 'sub.cluster'
endo <- RenameIdents(
  endo,
  "1_4" = "Trans",
  "11" = "Trans",
  "12" = "Trans",
  "15" = "Trans",
  "3" = "Trans",
  "6" = "Trans",
  "7" = "Trans",
  "8_16" = "Trans",
  "fcn1-0_14" = "FCN1+",
  "fcn1-13" = "FCN1+",
  "spp1-10" = "SPP1+MAM-",
  "spp1-rnase-2_9_0" = "Homeostatic",
  "spp1-rnase-2_9_1" = "SPP1+MAM-",
  "spp1-rnase-5" = "SPP1+MAM-"
)


endo <-
  SetIdent(endo, cells = Cells(subset(SPP1seu, tissue == "endo" &
                                        mam == "MAM")), value = "SPP1+MAM+")
endo$newclust <- Idents(endo)
endo$tissue <- 'endo'

table(endo$newclust)

# SPP1+MAM+       Trans       FCN1+   SPP1+MAM- Homeostatic
# 147        3348        1056         851         628


#####LUNG###########################
lung <- RenameIdents(
  lung,
  "0" = "Trans",
  "1_15" = "Trans",
  "10" = "Trans",
  "11" = "Trans",
  "14" = "Trans",
  "16" = "Trans",
  "3_12" = "Trans",
  "4_8_13" = "Homeostatic",
  "6" = "Trans",
  "9" = "Trans",
  "fcn1-7" = "FCN1+",
  "rnase-2" = "Homeostatic",
  "spp1-rnase-5" = "SPP1+MAM-"
)

lung <-
  SetIdent(lung, cells = Cells(subset(SPP1seu, tissue == "lung" &
                                        mam == "MAM")), value = "SPP1+MAM+")
lung$newclust <- Idents(lung)
lung$tissue <- 'lung'

table(lung$newclust)

# SPP1+MAM+       Trans Homeostatic       FCN1+   SPP1+MAM-
#   3189      107521       48446       11645       13496


#####HEART#################
heart <- RenameIdents(
  heart,
  "fcn1-6" = "FCN1+",
  "fcn1-4_10" = "FCN1+",
  "rnase-1" = "Homeostatic",
  "spp1-2" = "SPP1+MAM-",
  "spp1-rnase-9" = "SPP1+MAM-",
  "0" = "Trans",
  "3" = "Trans",
  "5" = "Trans",
  "7" = "Trans",
  "8" = "Trans"
)

heart <-
  SetIdent(heart, cells = Cells(subset(SPP1seu, tissue == "heart" &
                                         mam == "MAM")), value = "SPP1+MAM+")
heart$newclust <- Idents(heart)
heart$tissue <- 'heart'

table(heart$newclust)

# SPP1+MAM+       FCN1+ Homeostatic   SPP1+MAM-       Trans
# 364        4814        4790        3229       11264


#####LIVER############################
liver <- RenameIdents(
  liver,
  "fcn1-2" = "FCN1+",
  "rnase-1" = "Homeostatic",
  "spp1-7" = "SPP1+MAM-",
  "0" = "Homeostatic", #MACRO
  "3_5" = "Trans",
  "4" = "Trans",
  "6" = "Trans"
)

liver <-
  SetIdent(liver, cells = Cells(subset(SPP1seu, tissue == "liver" &
                                         mam == "MAM")), value = "SPP1+MAM+")
liver$newclust <- Idents(liver)
liver$tissue <- 'liver'

table(liver$newclust)

# SPP1+MAM+       FCN1+ Homeostatic   SPP1+MAM-       Trans
# 29        1680        4683         330        3679


# ### Keep only SPP1 and Trans populations
# skin <- DietSeurat(skin[ , skin$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                    assays = c("RNA","SCT","integrated"))
# kidney <- DietSeurat(kidney[ , kidney$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                      assays = c("RNA","SCT","integrated"))
# endo <- DietSeurat(endo[ , endo$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                    assays = c("RNA","SCT","integrated"))
# lung <- DietSeurat(lung[ , lung$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                    assays = c("RNA","SCT","integrated"))
# heart <- DietSeurat(heart[ , heart$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                     assays = c("RNA","SCT","integrated"))
# liver <- DietSeurat(liver[ , liver$newclust %in% c('SPP1+MAM+', 'SPP1+MAM-', 'Trans')], 
#                     assays = c("RNA","SCT","integrated"))

# obj.ls <- list(skin, kidney, endo, lung, heart, liver)
# SPP1Trans <- purrr::reduce(obj.ls, merge)

# downsample lung object
set.seed(42)
idx = sample.int(ncol(lung), size = 25000)
lung_ds <- lung[,idx]
lung_ds

table(lung_ds$newclust)
# SPP1+MAM+       Trans Homeostatic       FCN1+   SPP1+MAM-
#   404       14682        6476        1614        1824


# Combine all tissue into one object
obj.ls <- list(skin, kidney, endo, lung_ds, heart, liver)
all <- purrr::reduce(obj.ls, merge)

saveRDS(all, './data/all.rds')