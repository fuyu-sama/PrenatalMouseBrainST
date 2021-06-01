#! /usr/bin/env Rscript

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||  :  |||// \
#                / _||||| -:- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
#         \  \ `_.   \_ __\ /__ _/   .-` /  /
#     =====`-.____`.___ \_____/___.-`___.-'=====
#                       `=---='
#
#
#     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#               佛祖保佑         永无BUG
#  Codes are far away from bugs with Buddha's bless
#

# %% environment config
HOME <- Sys.getenv("HOME")
renv::activate(paste0(HOME, "/workspace/mouse-brain-full"))
library(Seurat)

arg <- commandArgs(trailingOnly = TRUE)
idx <- arg[1]
sessionInfo()

# %% read data and create Seurat object
count_df <- read.csv(
    paste0(
        HOME, "/workspace/mouse-brain-full/spaceranger/", idx,
        "/outs/filtered_feature_bc_matrix/", idx, ".csv"
    ),
    check.names = FALSE, row.names = 1
)
seurat_obj <- CreateSeuratObject(count_df, min.cells = 3)

# %% SCTransform
seurat_obj <- SCTransform(
    seurat_obj,
    return.only.var.genes = FALSE, verbose = FALSE
)

# %% dimension redunction and cluster
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
write.csv(
    Idents(seurat_obj),
    paste0(HOME, "/mouse-brain-full/SCT/SNN/",idx, ".csv")
)
