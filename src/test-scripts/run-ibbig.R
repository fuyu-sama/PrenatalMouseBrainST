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
library(iBBiG)

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

# %% run
res <- iBBiG(as.matrix(read_df), nModules = 15)
jpeg(
    paste0(WORKDIR, "results/3.jpg"),
    width = 1000, height = 1000
)
plot(res)
dev.off()
