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
library(SC3)
library(Seurat)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
set.seed(42)
sessionInfo()

# %% read data and create Seurat object
count_path <- paste0(
    HOME, "/workspace/mouse-brain-full/spaceranger/", idx,
    "/outs/filtered_feature_bc_matrix/", idx, ".csv"
)
save_path <- paste0(
    HOME, "/workspace/mouse-brain-full/SCT/SC3/pattern/",
    idx, "_SC3.csv"
)

#genes <- rownames(read.csv(paste0(HOME, "/workspace/mouse-brain-full/logcpm/scale_df/", idx, "-logcpm-inter.csv"),
#check.names = FALSE, row.names = 1))
count_df <- read.csv(count_path, check.names = FALSE, row.names = 1)
#count_df <- count_df[genes, ]
seurat_obj <- CreateSeuratObject(count_df)

# %% SCTransform
seurat_obj <- SCTransform(seurat_obj, return.only.var.genes = FALSE, verbose = FALSE)

# %% create SingleCellExperiment object and SC3
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(seurat_obj[["SCT"]]@counts),
        logcounts = as.matrix(seurat_obj[["SCT"]]@scale.data)
    )
)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(
    sce,
    ks = 5:28,
    biology = FALSE,
    gene_filter = FALSE,
    n_cores = 5
)

col_data <- colData(sce)
col_data <- data.frame(col_data[, grep("sc3_", colnames(col_data))])
write.csv(col_data, save_path)

