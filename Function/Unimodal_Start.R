library(Matrix)
library(Seurat)
library(cluster)
library(fitdistrplus)
library(data.table)
library(dplyr)
library(ggrepel)

source("/data/project/scitoseq/script/Final/Function/Make_hashtag_obj.R")
source("/data/project/scitoseq/script/Final/Function/Make_resolved_obj.R")

raw <- Read10X_h5('/data/project/scitoseq/script/Final/Sample/200K_28AB_10Pool/filtered_feature_bc_matrix.h5')
abNum <- 28
batchNum <- 10

hashtaglist <- make_hashtag_obj(raw = raw, abNum = abNum, batchNum = batchNum)
result <- make_resolved_obj(ab = hashtaglist$ab, tx = hashtaglist$tx, discrete = hashtaglist$discrete,
                                  ab_hashtag_obj = hashtaglist$ab_hashtag_obj,
                                  abNum = abNum, batchNum = batchNum)
selected_ab_cell_matrix <- after_processing(result$selected_ab_cell_matrix)
meta_ab_cell_matrix <- result$meta_ab_cell_matrix
ptn_names <- read.table("/data/project/scitoseq/script/Final/Sample/200K_28AB_10Pool/ptn_names.csv")

source("/data/project/scitoseq/script/Final/Function/python_uni.R")