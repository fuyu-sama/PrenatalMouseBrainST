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
WORKDIR <- paste0(Sys.getenv("HOME"), "/workspace/mouse-brain-full/")
renv::activate(WORKDIR)
library(SC3)
library(Seurat)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
set.seed(42)
sessionInfo()

# %% read data and create Seurat object
count_path <- paste0(WORKDIR, "scale_df/raw_count/", idx, "-raw.csv")
save_path <- paste0(WORKDIR, "scale_df/SCTransform/", idx, "-SCTransform.csv")

count_df <- read.csv(count_path, check.names = FALSE, row.names = 1)
seurat_obj <- CreateSeuratObject(count_df)

# %% SCTransform
seurat_obj <- SCTransform(
    seurat_obj,
    return.only.var.genes = FALSE, verbose = FALSE
)
write.csv(seurat_obj[["SCT"]]@scale.data, save_path)
