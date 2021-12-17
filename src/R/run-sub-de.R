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
library(Seurat)

options(future.globals.maxSize = 1024 * 1024^2)

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
scale_method <- args[1]
cluster_method <- args[2]
idx <- args[3]
region <- args[4]
ncs <- args[5]

# %% read data and create Seurat object
read_df <- read.csv(
    paste0(
        WORKDIR, "results/cluster/", scale_method, "-", cluster_method,
        "/sub-cluster/tables/", idx, "-", region, "-count.csv"
    ),
    check.names = F, row.names = 1
)
cluster_df <- read.csv(
    paste0(
        WORKDIR, "results/cluster/", scale_method, "-", cluster_method,
        "/sub-cluster/tables/", idx, "-", region, "-cluster.csv"
    ),
    check.names = F, row.names = 1
)
read_df <- read_df[, rownames(cluster_df)]
if (!all(rownames(cluster_df) == colnames(read_df))) {
    stop("AssertionError")
}

seurat_obj <- CreateSeuratObject(read_df)
seurat_obj <- SetIdent(seurat_obj, value = cluster_df[, ncs])

# %% DGE
de <- FindAllMarkers(seurat_obj, min.pct = 0.1)
write.csv(
    de,
    paste0(
        WORKDIR, "results/cluster/", scale_method, "-", cluster_method,
        "/sub-cluster/tables/", idx, "-", region, "-", ncs, "-de.csv"
    ))
