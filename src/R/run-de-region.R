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
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 1024 * 1024^2)
sessionInfo()

scale_method = "combat"

# %% read data and create Seurat object
read_df <- read.csv(
    paste0(
        WORKDIR,
        "Data/scale_df/", scale_method, "/full-", scale_method, "-inter.csv"),
    check.names = F, row.names = 1
)
cluster_df <- read.csv(
    paste0(WORKDIR, "results/cluster/", scale_method, "-SC3/pattern/full-SC3.csv"),
    check.names = F, row.names = 1
)
read_df <- read_df[, rownames(cluster_df)]
if (!all(rownames(cluster_df) == colnames(read_df))) {
    stop("AssertionError")
}

seurat_obj <- CreateSeuratObject(read_df)
seurat_obj <- SetIdent(seurat_obj, value = cluster_df[, 1])

regions <- jsonlite::read_json(
    paste0(WORKDIR, "results/cluster/", scale_method, "-SC3/regions.json"),
    simplifyVector = TRUE
)$regions

# %% DGE
de_1va_list <- list()
for (region in names(regions)) {
    print(region)
    de_1va_list[[region]] <- FindMarkers(
        seurat_obj,
        ident.1 = regions[[region]],
        min.pct = 0.1,
        verbose = TRUE
    )
    write.csv(
        de_1va_list[[region]],
        paste0(
            WORKDIR,
            "results/DE/", scale_method, "/region-specific/DE-", region, ".csv")
    )
}
