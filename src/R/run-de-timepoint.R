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

timepoints <- list(
    E135 = c("E135A", "E135B"),
    E155 = c("E155A", "E155B"),
    E165 = c("E165A", "E165B"),
    E175 = c("E175A1", "E175A2", "E175B"),
    P0 = c("P0A1", "P0A2")
)

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
scale_method <- args[1]
cluster_method <- args[2]

# %% read data and create Seurat object
read_df <- read.csv(
    paste0(
        WORKDIR,
        "Data/scale_df/", scale_method, "/full-", scale_method, ".csv"),
    check.names = F, row.names = 1
)
cluster_df <- read.csv(
    paste0(
        WORKDIR,
        "results/cluster/", scale_method, "-", cluster_method,
        "/pattern/full-", cluster_method, ".csv"
        ),
    check.names = F, row.names = 1
)
read_df <- read_df[, rownames(cluster_df)]
if (!all(rownames(cluster_df) == colnames(read_df))) {
    stop("AssertionError")
}

seurat_obj <- CreateSeuratObject(read_df)
seurat_obj <- SetIdent(seurat_obj, value = cluster_df[, 1])

regions <- jsonlite::read_json(
    paste0(WORKDIR, "results/cluster/", scale_method, "-", cluster_method, "/regions.json"),
    simplifyVector = TRUE
)$regions

# %% DGE
de_list <- list()
for (region in names(regions)) {
    for (timepoint in names(timepoints)) {
        g <- paste(region, timepoint, sep = "-")
        print(g)
        ident_1 <- c()
        ident_2 <- c()
        for (cluster in regions[[region]]) {
            s <- strsplit(cluster, "_")[[1]][1]
            if (s %in% timepoints[[timepoint]]) {
                ident_1 <- c(ident_1, cluster)
            } else {
                ident_2 <- c(ident_2, cluster)
            }
        }
        de_list[[g]] <- FindMarkers(
            seurat_obj,
            ident.1 = ident_1,
            ident.2 = ident_2,
            min.pct = 0.1,
            verbose = TRUE
        )
        write.csv(
            de_list[[g]],
            paste0(
                WORKDIR,
                "results/DE/", scale_method, "-", cluster_method, "/timepoint-specific/DE-",
                g, ".csv"
            )
        )
    }
}
