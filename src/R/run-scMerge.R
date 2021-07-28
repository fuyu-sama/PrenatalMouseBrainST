#! /usr/bin/env Rscript

#
#                       _oo0oo_
#                      o8888888o
#                      88" . "88
#                      (| -_- |)
#                      0\  =  /0
#                    ___/`---'\___
#                  .' \\|     |// '.
#                 / \\|||   =  |||// \
#                / _||||| - =- |||||- \
#               |   | \\\  -  /// |   |
#               | \_|  ''\---/''  |_/ |
#               \  .-\__  '-'  __/-. /
#             ___'. .'  /--.--\  `. .'___
#          ."" '<  `.___\_<|>_/___.' >' "".
#         | |  =  `- \`.;`\ _ /`;.`/ - `  = | |
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
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scMerge)

sessionInfo()

# %% read data
count_df <- read.csv(
    paste0(WORKDIR, "Data/scale_df/raw_count/full-raw-inter.csv"),
    check.names = FALSE,
    row.names = 1
)
logcpm_df <- read.csv(
    paste0(WORKDIR, "Data/scale_df/logcpm/full-logcpm-inter.csv"),
    check.names = FALSE,
    row.names = 1
)
select_genes <- intersect(rownames(count_df), rownames(logcpm_df))
count_df <- count_df[select_genes, ]
logcpm_df <- logcpm_df[select_genes, ]

batch_df <- data.frame(stringsAsFactors = TRUE)
for (i in colnames(count_df)) {
    batch_df <- rbind(batch_df, strsplit(i, "_")[[1]][1])
}
rownames(batch_df) <- colnames(count_df)
colnames(batch_df) <- "batch"
batch_df$batch <- as.factor(batch_df$batch)

sce <- SingleCellExperiment(
    list(counts = as.matrix(count_df), logcounts = as.matrix(logcpm_df)),
    colData = batch_df
)

data("segList", package = "scMerge")

# %% unsupervised scMerge
scMerge_unsupervised_all <- scMerge(
    sce_combine = sce,
    ctl = segList$mouse$mouse_scSEG,
    kmeansK = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10),
    assay_name = "scMerge_unsupervised_all",
    replicate_prop = 1
)
write.csv(
    assay(scMerge_unsupervised_all, "scMerge_unsupervised_all"),
    paste0(WORKDIR, "Data/scale_df/scMerge/full-scMerge-inter.csv")
)
