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

# %% read data and create Seurat object
read_df <- read.csv(
    paste0(WORKDIR, "Data/scale_df/logcpm/full-logcpm-inter.csv"),
    check.names = F, row.names = 1
)
cluster_df <- read.csv(
    paste0(WORKDIR, "results/cluster/SCT-SC3/pattern/full-SC3.csv"),
    check.names = F, row.names = 1
)
read_df <- read_df[, rownames(cluster_df)]
if (!all(rownames(cluster_df) == colnames(read_df))) {
    stop("AssertionError")
}

seurat_obj <- CreateSeuratObject(read_df)
seurat_obj <- SetIdent(seurat_obj, value = cluster_df[, 1])

# %% DGE
regions <- list(
    cortex = c(
        "E135A_1", "E135B_3", "E135B_9",
        "E155A_4", "E155B_5", "E155B_6",
        "E165A_1", "E165B_4",
        "E175A1_8", "E175A2_5", "E175A2_6", "E175B_6",
        "P0A1_5", "P0A2_1", "P0B_3"
        ),
    thalamus = c(
        "E135A_2", "E135B_6", "E135B_7",
        "E155A_5", "E155A_7", "E155A_8", "E155B_7",
        "E165A_3", "E165A_5", "E165B_3",
        "E175A1_2", "E175A1_7", "E175A2_10", "E175A2_11", "E175A2_12",
        "E175B_2", "E175B_4", "E175B_5",
        "P0A1_6", "P0A1_7", "P0A1_12", "P0A2_6",
        "P0B_1", "P0B_7", "P0B_8"
        ),
    hypothalamus = c(
        "E135A_3", "E135A_6", "E135B_5", "E135B_12",
        "E155A_6", "E155B_1",
        "E165A_6", "E165B_1",
        # "E175A1_1", "E175A2_1",
        "E175B_3",
        "P0A1_1", "P0A2_4", "P0B_12"
        ),
    olfactory = c(
        "E135A_5", "E135A_7", "E135B_4",
        "E155A_11", "E155B_2", "E155B_11",
        "E165A_10", "E165A_11", "E165B_7", "E165B_10",
        "E175A1_5", "E175A2_7", "E175A2_9", "E175B_1", "E175B_9",
        "P0A1_10", "P0A2_7", "P0B_5", "P0B_6"
    )
)

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
        paste0(WORKDIR, "results/DE/region-specific/DE-", region, ".csv")
    )
}
