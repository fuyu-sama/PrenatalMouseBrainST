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
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)
set.seed(16)
sessionInfo()

idx <- args[1]
scale_method <- args[2]

# %% read data and create sce
count_path <- paste0(WORKDIR, "Data/scale_df/raw/", idx, "-raw.csv")
sct_path <- paste0(
    WORKDIR, "Data/scale_df/", scale_method, "/", idx, "-", scale_method, ".csv")
save_path <- paste0(WORKDIR, "results/cluster/", scale_method, "-sc3/pattern/", idx, "-sc3.csv")
count_df <- read.csv(count_path, check.names = FALSE, row.names = 1)
sct_df <- read.csv(sct_path, check.names = FALSE, row.names = 1)
count_df <- count_df[rownames(sct_df), ]
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(count_df),
        logcounts = as.matrix(sct_df)
    )
)
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]

# %% sc3
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
