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
gene_list <- args[3]
suffix <- args[4]

if (is.na(idx)) idx <- "E165A"
if (is.na(scale_method)) scale_method <- "combat"

# %% read data and create sce
count_path <- paste0(WORKDIR, "Data/scale_df/raw/", idx, "-raw.csv")
scale_path <- paste0(
    WORKDIR, "Data/scale_df/", scale_method, "/", idx, "-", scale_method, ".csv")

count_df <- read.csv(count_path, check.names = FALSE, row.names = 1)
scale_df <- read.csv(scale_path, check.names = FALSE, row.names = 1)
count_df <- count_df[rownames(scale_df), ]

if (! is.na(gene_list)) {
    conn <- file(gene_list, open = "rt")
    gene_list <- readLines(conn)
    close(conn)
    count_df <- count_df[gene_list, ]
    scale_df <- scale_df[gene_list, ]
    scale_method <- paste(scale_method, suffix, sep = "-")
}
sce <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(count_df),
        logcounts = as.matrix(scale_df)
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
save_path <- paste0(
    WORKDIR, "results/cluster/", scale_method, "-sc3/pattern/", idx, "-sc3.csv")
write.csv(col_data, save_path)
