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
library(dplyr)
library(SC3)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
set.seed(42)
sessionInfo()

idx <- args[1]
scale_method <- args[2]
region <- args[3]
if (is.na(idx) | is.na(scale_method) | is.na(region)) {
    idx <- "E165A"
    scale_method <- "combat"
    region <- "cortex"
}

# %% read data
count_path <- paste0(WORKDIR, "Data/scale_df/raw/", idx, "-raw.csv")
sct_path <- paste0(
    WORKDIR, "Data/scale_df/", scale_method, "/", idx, "-", scale_method, ".csv")
count_df <- read.csv(count_path, check.names = FALSE, row.names = 1)
scale_df <- read.csv(sct_path, check.names = FALSE, row.names = 1)

# %% read cluster results
cluster_df <- read.csv(
    paste0(
        WORKDIR,
        "results/cluster/", scale_method, "-sc3",
        "/pattern/full-sc3.csv"
        ),
    check.names = F, row.names = 1
)
cluster_df <- as.data.frame(cluster_df[colnames(count_df), ])
rownames(cluster_df) <- colnames(count_df)
colnames(cluster_df) <- "sc3_clusters"

regions <- jsonlite::read_json(
    paste0(WORKDIR, "results/cluster/", scale_method, "-sc3/regions.json"),
    simplifyVector = TRUE
    )$regions

# %% subclustering
region_label <- c()
for (i in regions[["cortex"]]) {
    if (grepl(idx, i)) region_label <- c(region_label, i)
}

cluster_sub <- filter(cluster_df, sc3_clusters %in% region_label)
scale_sub <- scale_df[, rownames(cluster_sub)]
count_sub <- count_df[, rownames(cluster_sub)]

sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(count_sub),
        logcounts = as.matrix(scale_sub)
    )
)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

sce <- sc3(
    sce,
    ks = 2:9,
    biology = FALSE,
    gene_filter = FALSE,
    n_cores = 5
)

col_data <- as.data.frame(colData(sce))
write.csv(
    col_data,
    paste0(
        WORKDIR, "results/cluster/", scale_method,
        "-sc3/sub-clusters/pattern", idx, "-", region, "-sc3.csv"
    )
)
