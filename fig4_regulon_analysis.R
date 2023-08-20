setwd('C:/Users/e0205142/OneDrive - National University of Singapore/OtherPj/Kevin/fig4/')
library(Seurat)
library(SCENIC)
library(SCopeLoomR)
library(tidyverse)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(data.table)
library(rvg)
library(officer)
library(ggrastr)
library(ggpubr)

SPP1seu <- readRDS('./data/SPP1seu.rds')
SPP1seu

SPP1.loom <- open_loom('./final_output/scenic/SPP1seu_SCENIC_res.loom',mode = 'r')

# Regulon activity -------------
## 1. retrieve regulon activity matrix (AUC) and thresholds from loom file
AUC_SPP1 <- get_regulons_AUC(SPP1.loom, column.attr.name='RegulonsAUC' )
Thresholds_SPP1 <- get_regulon_thresholds(SPP1.loom)
regulonsAucThresholds_SPP1 <- setNames(as.numeric(names(Thresholds_SPP1)), Thresholds_SPP1)

## 2. Binarize regulon activity matrix
regulons_activity_SPP1 <- AUC_SPP1@assays@data$AUC %>% as.data.frame() %>%
  mutate(across(everything(), ~ case_when(
    .x > regulonsAucThresholds_SPP1 ~ 1,
    .x <= regulonsAucThresholds_SPP1 ~ 0
  )))

dim(regulons_activity_SPP1)
# [1]   238 24105


## We only keep regulons activated in >= 10% cells in at least 4 tissues 
## as described in Methods session.

# transpose regulons_activity_SPP1
regulons_activity_SPP1_t <- transpose(regulons_activity_SPP1)
rownames(regulons_activity_SPP1_t) <- colnames(regulons_activity_SPP1)
colnames(regulons_activity_SPP1_t) <- rownames(regulons_activity_SPP1)

# kept_regulons_tmp <- regulons_activity_SPP1_t %>% mutate(tissue = SPP1seu$tissue) %>%
#   pivot_longer(cols = -tissue, names_to = 'TF', values_to = 'activity') %>% group_by(tissue, TF) %>%
#   summarise(act_ratio = mean(activity), .groups = 'drop') %>% group_by(TF) %>% 
#   filter(sum(act_ratio>=0.1) >= 4) %>% pull(TF) %>% unique()

kept_regulons <- regulons_activity_SPP1_t %>% mutate(tissue = SPP1seu$tissue) %>%
  group_by(tissue) %>% summarise(across(everything(), ~ mean(. > 0))) %>% 
  select(where(function(x) sum(x >= 0.1) >= 4)) %>% select(-tissue) %>% 
  colnames()

length(kept_regulons)
# [1] 173

# Fig4A-------------------------
## DOR calculation -------------------------

## merge metadata and regulons activity to the same dataframe
DOR_df_input <- regulons_activity_SPP1_t %>% 
  select(kept_regulons) %>%
  mutate(tissue = SPP1seu$tissue,
         mam = SPP1seu$mam)

## define DOR calculation function
cal_DOR <- function(mam_state, regulon_activity){
  regulon_activity <- factor(regulon_activity, levels = c(0,1))
  conf_matrix <- table(mam_state, regulon_activity)
  #print(conf_matrix)
  TP = conf_matrix['MAM', '1']
  FN = conf_matrix['MAM', '0']
  TN = conf_matrix['nonMAM', '0']
  FP = conf_matrix['nonMAM', '1']
  DOR = ((TP+.05)*(TN+.05))/((FP+.05)*(FN+.05))
  return(DOR)
}

## split data by tissue and calculate DOR for each tissue
tissue_gp <- DOR_df_input %>%
  dplyr::group_by(tissue) %>% group_split()


regulons_DOR_list <- list()
for (df in tissue_gp) {
  regulons_DOR <- sapply(df[,1:(ncol(df)-2)], cal_DOR, mam_state=df$mam)
  regulons_DOR_list[[unique(df$tissue)]] <- regulons_DOR
}

regulons_DOR_df <- data.frame(regulons_DOR_list)
regulons_DOR_df

# check distribution of DOR in each tissue
ggData <- regulons_DOR_df %>% pivot_longer(cols = everything(), 
                                           names_to = 'tissue',
                                           values_to = 'DOR')
head(ggData)
ggplot(ggData, aes(DOR))+
  geom_histogram(bins = 50)+
  facet_wrap(~tissue, scales = "free")

###############################-----------------###
### Mean value excluding outliers
library(outliers) #containing function outlier

regulons_DOR_mean <- regulons_DOR_df %>% 
  rownames_to_column(var='TF') %>% 
  pivot_longer(cols = -TF, names_to = 'tissue') %>%
  group_by(TF) %>%
  filter(!value %in% c(outlier(value))) %>%
  summarise(mean_DOR = mean(value, na.rm = TRUE)) %>%
  arrange(desc(mean_DOR))
#############################-------------------###

### Median of DOR
regulons_DOR_median <- rowMedians(as.matrix(regulons_DOR_df))

regulons_DOR_median_df <- data.frame('DOR' = regulons_DOR_median,
                                     'log10_DOR' = log10(regulons_DOR_median),
                                     row.names = rownames(regulons_DOR_df))

regulons_DOR_median_df %<>% arrange(desc(DOR)) %>% mutate(Rank=1:nrow(regulons_DOR_median_df))

regulons_DOR_median_df$TF_names <- gsub("\\([+-]\\)", "", rownames(regulons_DOR_median_df))
head(regulons_DOR_median_df)

## Plot------------------
ggData_4a <- regulons_DOR_median_df %>% rownames_to_column(var='TF')
head(ggData_4a)
p4a <- ggplot(data = ggData_4a, aes(x = Rank, y = log10_DOR, label = TF)) + 
  geom_point(size=0.7) + 
  theme_classic(base_size = 20) +
#  geom_hline(yintercept = log10(10), color = 'red',linetype='dashed') + 
  ylab('DOR for regulon activation\n in SPP1+MAM+ (log10)') + xlab('Rank') +
  ggrepel::geom_text_repel(data = subset(ggData_4a, Rank<=20), 
                           aes(x = Rank, y = log10_DOR), 
                           size = 4, seed = 42,
                           force        = 0.032,
                           force_pull   = 0,
                           nudge_x      = 150,
                           vjust = 0,
                           direction    = "y")+
  # ggrepel::geom_text_repel(data = subset(ggData_4a, DOR<=0.25), 
  #                          aes(x = Rank, y = log10_DOR), 
  #                          size = 5, seed = 42,
  #                          force        = 0.028,
  #                          force_pull   = 0,
  #                          nudge_x      = -50,
  #                          vjust = 0,
  #                          direction    = "y")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=22))

# ori <- c('NFATC1', 'HIVEP3', 'MAF', 'NR1H3', 'ETV5', 'THRA')
# lapply(ori,grep,x=ggData$TF_names)
p4a

ggsave(p4a, filename='./final_output/fig/fig4a.pdf', width = 7, height = 8)

# ori <- c('NFATC1', 'HIVEP3', 'MAF', 'NR1H3', 'ETV5', 'THRA')
# p4a_ori <- ggplot(data = ggData, aes(x = Rank, y = log10_DOR, label = TF_names)) + 
#   geom_point(size=0.7) + 
#   theme_classic(base_size = 25, base_line_size = 0.1) +
#   geom_hline(yintercept = 1, color = 'red') + ylab('DOR (log10)') + xlab('Rank of TF') +
#   ggrepel::geom_text_repel(data = subset(ggData, DOR>=10), aes(x = Rank, y = log10_DOR), size = 4, hjust=1,vjust=1) +
#   geom_point(data = subset(ggData, TF_names%in%ori), color='red', size=2) + 
#   ggrepel::geom_text_repel(data = subset(ggData, TF_names%in%ori), aes(x = Rank, y = log10_DOR), size = 4, hjust=1,vjust=1.5, color='red')
# p4a_ori
# ggsave(p4a_ori, filename='./final_output/fig4a_ori.pdf', width = 5, height = 5)

# Fig4B---------------------------------
## MAM signature loading-----------
all_marker <- read.delim('./data/mam_DEsumm.txt')
rownames(all_marker) <- all_marker$gene
head(all_marker)

marker <- all_marker$gene[all_marker$avgLFC >0.5 & all_marker$zDiffProp >1.96]
marker_d <- all_marker$gene[all_marker$avgLFC < -0.5 & all_marker$zDiffProp < -1.96]

## Regulons Specificity calculation-----------------
## regulon specificity to mam signature

head(regulons_DOR_median_df) #(row=TF, column=target genes)

regulons_SPP1 <- get_regulons(SPP1.loom, column.attr.name='Regulons')
head(regulons_SPP1[1:10,1:10])

dim(regulons_SPP1)
# [1]  238 3000
dim(regulons_DOR_median_df)
# [1] 173  4

### target genes-------------
regulons_genes_SPP1 <- apply(regulons_SPP1, 1, function(x){names(x)[x==1]})
regulons_genes_SPP1 <- regulons_genes_SPP1[rownames(regulons_DOR_median_df)]
length(regulons_genes_SPP1)

int_u <- c()
int_d <- c()
reg_l <- c()
reg_d <- c()
reg_len <- c()

for(reg in regulons_genes_SPP1) {
  q <- length( intersect( marker, reg)) - 1
  m <- length(marker)
  n <- nrow(all_marker) - m
  reg_l <- c(reg_l, phyper(q, m, n, length(reg), lower.tail = FALSE, log.p = FALSE))
  
  q <- length( intersect( marker_d, reg)) - 1
  m <- length(marker)
  n <- nrow(all_marker) - m
  reg_d <- c(reg_d, phyper(q, m, n, length(reg), lower.tail = FALSE, log.p = FALSE))
  
  reg_len <- c(reg_len, length(reg))
  int_u  <- c(int_u , sum(all_marker[intersect( marker, reg), 'avgLFC']))
  
  int_d  <- c(int_d , sum(all_marker[intersect( marker_d, reg), 'avgLFC'] ))
}

head(int_u) # sum of logFC in up-regulated markers (positive value)
head(int_d) # sum of logFC in down-regulated markers (negative value)
head(reg_len) # p_value of hypergeometric test for up markers
head(reg_l) # p_value of hypergeometric test for down markers

names(reg_len) <- names(regulons_genes_SPP1)
names(reg_l) <- names(regulons_genes_SPP1)
names(int_u) <- names(regulons_genes_SPP1)
names(int_d) <- names(regulons_genes_SPP1)
names(reg_d) <- names(regulons_genes_SPP1)
reg_l <- log10(p.adjust(reg_l)) * -1
reg_d <- log10(p.adjust(reg_d)) * -1

final <- (abs(reg_l - reg_d) / reg_len * (int_u + int_d)) #[ rowSums(prop_ALL > .1) >= 6]

# Expression of regulons (module score, x-axis in fig4b)
### ModuleScore in potato-----------
DefaultAssay(SPP1seu) <- 'SCT'
#DefaultAssay(SPP1seu) <- 'integrated'
## combine selected regulon and mam signature into one list
all_reg <- regulons_genes_SPP1
all_reg[['SPP1+MAM+']] <- marker
saveRDS(all_reg, './final_output/rds/all_regulon_mamSig.rds')

SPP1seu <- AddModuleScore(SPP1seu, features = all_reg, 
                          name = "_modulescore")
names(SPP1seu@meta.data)[grep("_modulescore", names(SPP1seu@meta.data))] <- names(all_reg)
colnames(SPP1seu@meta.data)

# Extract regulon module score of each cell from metadata
regulons_expr <- SPP1seu@meta.data[, names(regulons_genes_SPP1)]
regulons_expr$mam <- SPP1seu$mam
head(regulons_expr)

# Calculate difference in regulon expression (SPP1+MAM+ minus SPP1+MAM-)
regulons_expr_df <- regulons_expr %>% 
  pivot_longer(-mam, names_to = 'row.names') %>% 
  group_by(mam, row.names) %>% 
  summarise(mean = mean(value)) %>% ungroup() %>% 
  pivot_wider(names_from = mam, values_from = mean) %>% 
  mutate(exprDiff = MAM-nonMAM) %>% as.data.frame()

regulons_DOR_median_df$row.names <- rownames(regulons_DOR_median_df)

head(regulons_DOR_median_df)
head(regulons_expr_df)

## all results in one table
regulons_DOR_median_df  <- left_join(regulons_DOR_median_df, regulons_expr_df, by = 'row.names')
regulons_DOR_median_df$Specificity <- final[regulons_DOR_median_df$row.names]
regulons_DOR_median_df$regulon_size <- reg_len[regulons_DOR_median_df$row.names]

head(regulons_DOR_median_df)

## plot------------------
ggData_4b <- regulons_DOR_median_df

## x-axis = difference in regulon expression (module score of target genes)
# p4b <- ggplot(ggData, aes(exprDiff, Specificity, label = TF_names)) +
#   geom_point(shape = 21, stroke = .1, aes(size = regulon_size), color = 'black') +
#   #geom_point(data = subset(ggData, DOR>=10), shape = 21, stroke = .1, aes(size = regulon_size), color = 'red')+
#   #geom_point(data = subset(ggData, DOR<10), shape = 21, stroke = .1, aes(size = regulon_size), color = 'black') + 
#   scale_size(range = c(2,20)) +
#   geom_text_repel(data = subset(ggData, Specificity>=2.5), 
#                   seed = 1, size = 5) +
#   theme_classic(base_size = 20, base_line_size = .5) + 
#   ylab('Regulon Specificity') + xlab('Difference in regulon expression') +
#   geom_vline(xintercept  = 0, linetype = "longdash", colour="#D3D3D3") +
#   geom_hline(yintercept  = 0, linetype = "longdash", colour="#D3D3D3")

p4b <- ggplot(ggData_4b, aes(exprDiff, Specificity, label = row.names)) +
  geom_vline(xintercept  = 0, linetype = "longdash", colour="#D3D3D3") +
  geom_hline(yintercept  = 0, linetype = "longdash", colour="#D3D3D3")+
  geom_point(size=2, color = 'black', alpha=0.2) +
  geom_point(data = subset(ggData_4b, Rank<=20&Specificity>0&exprDiff>0), size=2, color = 'red')+
  #geom_point(data = subset(ggData_4b, Specificity < -0.3), size=2, color = 'red') + 
  scale_size(range = c(2,20)) +
  geom_text_repel(data = subset(ggData_4b, Rank<=20&Specificity>0&exprDiff>0), 
                  seed = 1, size = 5)  +
  #geom_text_repel(data = subset(ggData_4b, Specificity < -0.3), 
  #                seed = 1, size = 5)+ # show bottom top 5
  theme_classic(base_size = 20) + 
  ylab('Regulon specificity to \nSPP1+MAM+ signature') + 
  xlab('Difference in regulon expression\n(SPP1+MAM+ - SPP1+MAM-)')+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=22))

p4b
ggsave(p4b, filename='./final_output/fig/fig4b.pdf', width = 10, height = 8)


# result dataframe---------
head(regulons_DOR_median_df)
write.csv(regulons_DOR_median_df, './final_output/table/fig4_regulons_scores.csv',
          row.names = F)

# Candidates------------------
selected_regulons_df <- regulons_DOR_median_df %>% filter(DOR>=10&Specificity>0&exprDiff>0)
selected_regulons_df

selected_r_names <- selected_regulons_df$TF_names
selected_r <- selected_regulons_df$row.names
select_r_list <- regulons_genes_SPP1[selected_r]
names(select_r_list) <- selected_r_names

selected_r_names
selected_r
select_r_list

save(selected_r_names, selected_r, select_r_list, 
     file = "./final_output/rds/selected_regulons.Rdata")

#### export regulon genes ###
# Load the openxlsx package
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()

# Add a worksheet to the workbook
addWorksheet(wb, "Sheet1")

# Write the vectors into separate columns
for (i in seq_along(regulons_genes_SPP1)) {
  col_num <- i  # Column number to write the vector
  col_name <- names(regulons_genes_SPP1)[i]  # Column name
  vector <- regulons_genes_SPP1[[i]]  # Vector to write
  
  # Write the column name
  writeData(wb, "Sheet1", col_name, startCol = col_num, startRow = 1)
  
  # Write the vector values below the column name
  writeData(wb, "Sheet1", vector, startCol = col_num, startRow = 2)
}

# Save the workbook as an Excel file
saveWorkbook(wb, "./final_output/table/SuppT9_Member genes of regulons.xlsx", overwrite = TRUE)

# AUCell on all cells-------------------
## refer to create_all_object.R to get all.rds
## run following code on HPC
# cd /data/petretto/home/e0205142/project/Kevin/fig4
# conda activate r4
# R
library(Seurat)
library(AUCell)

all <- readRDS('../data/all.rds')
all.ls <- SplitObject(all, split.by = 'tissue')

all_reg <- readRDS('./final_output/rds/all_regulon_mamSig.rds')

### WATCH OUT! use same assay as input to pySCENIC for AUC calculation
all_AUC_list <- list()
for (tissue in names(all.ls)){
  tissue_seu <- all.ls[[tissue]]
  expr_ranked <- AUCell_buildRankings(tissue_seu$integrated@data, plotStats=FALSE)
  sig_AUC_res <- AUCell_calcAUC(all_reg, expr_ranked)
  sig_AUC <- as.data.frame(t(sig_AUC_res@assays@data$AUC))
  # Add meta data inside
  sig_AUC$tissue <- tissue
  sig_AUC$clust <- tissue_seu$clust
  sig_AUC$newclust <- tissue_seu$newclust
  
  all_AUC_list[[tissue]] <- sig_AUC
  print(tissue)
}

all_AUC_df <- do.call("rbind", all_AUC_list)
saveRDS(all_AUC_list, './final_output/rds/all_AUC_list_inte_based.rds')
saveRDS(all_AUC_df, './final_output/rds/all_AUC_df_inte_based.rds')
#write.csv(all_AUC_df, './final_output/all_AUC_df_inte_based.csv')

# Fig4D--------------------------------------------
## plot regulon module score & SPP1+MAM+ signature module score
# color to use
DefaultAssay(SPP1seu) <- 'SCT'

# library(scales)

# test <- apply(SPP1seu@meta.data[,selected_r], 2, rescale)
# test <- as.data.frame(test)
# new_col <- selected_r_names
# colnames(test) <- new_col
# SPP1seu <- AddMetaData(SPP1seu, test)
# FeaturePlot(SPP1seu, features = selected_r_names, raster = F,raster.dpi = c(1024, 1024),
#             pt.size = 0.0001, combine = T, order = T)

FeaturePlot(SPP1seu, features = selected_r, raster = F,raster.dpi = c(1024, 1024),
            pt.size = 0.0001, combine = T, max.cutoff='q95', keep.scale = 'all')

color_min = min(SPP1seu@meta.data[,selected_r])
color_max = max(SPP1seu@meta.data[,selected_r])

p_expr <- FeaturePlot(SPP1seu, features = c(selected_r, 'SPP1+MAM+'), raster = F,raster.dpi = c(1024, 1024),
                      pt.size = 0.0001,combine = F,max.cutoff='q95')

#DotPlot(SPP1seu, features = c('JDP2', 'KLF3', 'CEBPD')) & scale_color_distiller(palette = 'RdYlBu')

p_expr_2 <- list()
i=1
for (p in p_expr){
  print(i)
  p_expr_2[[i]] <- p+scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))[c(2,3,4,5,6,7,8,9,10,11,11)],limits=c(color_min, color_max)) +
    theme(aspect.ratio = 1)
  #+labs(title = selected_r_names[i])
  i=i+1
}
patchwork::wrap_plots(p_expr_2, ncol=2)

p_umap_mam <- DimPlot(SPP1seu, group.by = 'mam', raster = F,raster.dpi = c(1024, 1024),
                      pt.size = 0.0001, shuffle = T) + theme(aspect.ratio = 1) +
  labs(title = 'Cell State') +
  scale_color_manual(labels = c("SPP1+MAM-", "SPP1+MAM+"), values = c("brown","grey"))
p_umap_mam

p_expr_2[[4]] <- p_umap_mam

# wrap potatoes together and remove axis
p4d <- patchwork::wrap_plots(p_expr_2, ncol=2) + plot_layout(guides = "collect")
p4d <- p4d&theme(line = element_blank(),
                 axis.title = element_blank(),
                 axis.text=element_blank(),
                 axis.ticks=element_blank()
)


p4d <- rasterise(p4d,dpi = 60)
p4d

ggsave(p4d,
       filename='./final_output/fig/fig4d.pdf', height = 12, width = 12)


# Fig4C-----------------
# load AUC for regulons and mam signatures in all cells
all_AUC_df <- readRDS('./final_output/rds/all_AUC_df_inte_based.rds')

sig_AUC_df <- all_AUC_df[!(all_AUC_df$newclust%in%c('FCN1+','Homeostatic')),c(selected_r, 'SPP1+MAM+', 'tissue', 'newclust')]
head(sig_AUC_df)
### DOWNSAMPLE  FOR VISUALIZATION ####
set.seed(8)

sig_AUC_df_down <- sig_AUC_df[sample(rownames(sig_AUC_df), size  = 10000, replace = F), ] 
#sig_AUC_df_down <- sig_AUC_df[sig_AUC_df$tissue=='lung',]

head(sig_AUC_df_down)
table(sig_AUC_df_down$tissue, sig_AUC_df_down$newclust)

ggData <- sig_AUC_df_down %>% 
  pivot_longer(cols = 1:length(selected_r), names_to = 'TF')
ggData$TF <- factor(ggData$TF, levels = selected_r)
head(ggData)

p4c <- ggplot(data = ggData, aes(x = `SPP1+MAM+`, y = value, color = newclust)) + 
  #ggrastr::geom_point_rast(alpha = .2, size = 1.2) +
  geom_point(alpha = .2, size = 1.2) +
  stat_smooth(data = subset(ggData, newclust %in% c('Trans', 'SPP1+MAM-')), 
              method="lm", se=.99995, fill=NA,  formula=y ~ poly(x, 1, raw=TRUE),colour="black")   + 
  stat_smooth(data = subset(ggData, newclust %in% c('SPP1+MAM-', 'SPP1+MAM+')),
              method="lm", se=.99995, fill=NA,  formula=y ~ poly(x, 1, raw=TRUE),colour="red") +
  facet_wrap(~TF, ncol = 2,scales='free') +
  theme(axis.line=element_line(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=22)) + 
  theme_classic(base_size = 20)  + 
  ylab('Regulon activity score (AUC)') + 
  xlab('SPP1+MAM+ activity score (AUC)')  + 
  scale_shape(solid = FALSE)+ 
  theme(aspect.ratio = 1, strip.background = element_blank(),
        panel.spacing.x=unit(2, "lines"),
        axis.text = element_text(size = 12))+
  guides(colour = guide_legend(override.aes = list(size=3,alpha = 1)))

p4c

ggsave(p4c, filename='./final_output/fig/fig4c.pdf', 
       width = 9, height = 9)

# Fig4E----------------------------
## plot AUC of selected regulons tissue by tissue
all_AUC_df <- readRDS('./final_output/rds/all_AUC_df_inte_based.rds')

fig4e.data.list <- list()

for (t in unique(all_AUC_df$tissue)){
  all_AUC_df.sub <- all_AUC_df %>% 
    filter(tissue == t) %>%
    filter(newclust %in% c('Homeostatic','Trans', 'SPP1+MAM-', 'SPP1+MAM+' ))
  
  median_home <- all_AUC_df.sub %>% select(newclust, selected_r) %>%
    filter(newclust=='Homeostatic') %>% select_if(is.numeric) %>%
    summarize(across(everything(), median))
  
  mscore_ratio <- sweep(all_AUC_df.sub[, selected_r], 2, as.numeric(median_home[1,]), '/')
  mscore_ratio$newclust <- all_AUC_df.sub$newclust

  ggData_sub <- mscore_ratio %>% pivot_longer(cols = -newclust, names_to = 'TF')
  ggData_sub$newclust <- factor(ggData_sub$newclust, levels = c('Homeostatic','Trans',  'SPP1+MAM-', 'SPP1+MAM+' ))
  ggData_sub$tissue <- t
  
  fig4e.data.list[[t]] <- ggData_sub
  
}

ggData.4e <- do.call("rbind", fig4e.data.list)
ggData.4e$TF <- fct_rev(factor(ggData.4e$TF, levels = selected_r))

p4e <- ggplot(ggData.4e, aes(x=TF, y=value, fill=newclust))+
  geom_boxplot(lwd=.2) + 
  #geom_hline(yintercept=10)+
  facet_grid(. ~ tissue, scales = 'free')+
  theme_classic(base_size = 20)+
  #coord_flip(ylim = c(-1,10))+
  coord_flip()+
  guides(fill = guide_legend(reverse = TRUE))+
  ylab('Fold change of regulon activity (AUC)\n in comparison with the median expression in homeostatic macrophages') + 
  xlab('Regulons')+
  theme(strip.background = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=22))+
  scale_fill_discrete(labels=c('Homeostatic','Transitional',  'SPP1+MAM-', 'SPP1+MAM+'))

p4e
ggsave(p4e, filename='./final_output/fig/fig4e_auc.pdf', height = 6, width=15)

## statistical test for fig4e
res <- compare_means(value~newclust, ggData.4e, group.by=c("TF", "tissue"), method = "wilcox.test")
res <- res %>% arrange(desc(p.adj))

write.csv(res, './final_output/table/fig4e_wilcoxon_test_result_AUCell.csv')


# Figure 4---------------------
pw <- ((p4a+p4b)/(p4c+p4d))/p4e
# pw <- pw&theme(axis.text=element_text(size=15),
#                axis.title=element_text(size=20))
#+ plot_annotation(tag_levels = 'A')
ggsave(pw,
       filename='./final_output/fig/fig4_r_output.pdf', height = 20, width=15)


