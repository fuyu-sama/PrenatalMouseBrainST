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
library(biclust)

options(future.globals.maxSize = 1024 * 1024^2)
sessionInfo()

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
scale_method <- args[2]
cluster_method <- args[3]
if (is.na(idx) | is.na(scale_method) | is.na(cluster_method)) {
    idx <- "E165A"
    scale_method <- "combat-zjq"
    cluster_method <- "sc3"
}

# %% read data
read_df <- read.csv(
    paste0(
        WORKDIR,
        "Data/scale_df/", scale_method, "-moran/",
        idx, "-", scale_method, "-moran.csv"),
    check.names = F, row.names = 1
)

# %% BCBimax
result <- biclust(
    as.matrix(read_df),
    method = BCBimax(),
    number = 10, minr = 50, minc = 20
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCBimax.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)

# %% BCCC
# H(I, J) = Sigma_(i in I, j in J)(a_i,j - a_i,J - a_I,j + a_IJ) ^ 2 / ||I||||J||
# delta: Maximum of accepted score; alpha: Scaling factor
result <- biclust(
    as.matrix(read_df),
    method = BCCC(),
    number = 10, delta = 1, alpha = 1
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCCC.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)

# %% BCPlaid
result <- biclust(
    as.matrix(read_df),
    method = BCPlaid(),
    verbose = FALSE
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCPlaid.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)

# %% BCQuest
result <- biclust(
    as.matrix(read_df),
    method = BCQuest()
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCQuest.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)

# %% BCSpectral
result <- biclust(
    as.matrix(read_df),
    method = BCSpectral()
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCSpectral.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)

# %% BCXmotifs
result <- biclust(
    as.matrix(read_df),
    method = BCXMotifs()
)

writeBiclusterResults(
    paste0(WORKDIR, "results/4/BCXmotifs.csv"), result,
    bicName = capture.output(str(result@Parameters$Call)), delimiter = ",",
    geneNames = rownames(read_df), arrayNames = colnames(read_df)
)
