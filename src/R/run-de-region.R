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
library(dplyr)
library(ggplot2)

options(future.globals.maxSize = 1024 * 1024^2)
sessionInfo()

# %% read data and create Seurat object
read_df <- read.csv(
    paste0(
        HOME, "/workspace/mouse-brain-full/",
        "/scale_df/logcpm/full-logcpm-inter.csv"
        ),
    check.names = F, row.names = 1
)
cluster_df <- read.csv(
    paste0(HOME, "/workspace/mouse-brain-full/results/SCT-SC3/full-12-SC3.csv"),
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
        "E135A_1", "E135B_11",
        "E155A_3", "E155B_9",
        "E165A_1", "E165B_10",
        "E175A1_4", "E175A2_2", "E175B_6",
        "P0A1_6", "P0A2_1", "P0B_1", "P0B_3"
        ), # 皮层
    hippocampus = c(
        "E135A_9", "E135B_2",
        "E155A_2", "E155B_8",
        "E165A_2", "E165B_9",
        "E175A1_6", "E175A2_5", "E175A2_6", "E175B_9",
        "P0A1_9", "P0A2_6", "P0B_4"
        ), # 海马
    thalamus = c(
        "E135A_2", "E135B_10",
        "E155A_12", "E155B_10",
        "E165A_3", "E165B_3",
        "E175A1_1", "E175A1_11", "E175A2_4", "E175A2_11", "E175B_5",
        "P0A1_8", "P0A2_5", "P0B_6", "P0B_7"
        ), # 丘脑
    hypothalamus = c(
        "E135A_6", "E135A_7", "E135B_3", "E135B_9",
        "E155A_11", "E155B_2",
        "E165A_4", "E165A_12", "E165B_1",
        "E175A1_9", "E175A2_10", "E175B_4",
        "P0A1_2", "P0A2_3", "P0B_12"
    ) # 下丘脑
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
        paste0(
            HOME, "/workspace/mouse-brain-full/",
            "logcpm/DE/region-specific/DE-10-", region, ".csv"
        )
    )
}

save.image(
    paste0(HOME, "/workspace/mouse-brain-full/",
        "logcpm/DE/region-specific/DE.RData"
    )
)
