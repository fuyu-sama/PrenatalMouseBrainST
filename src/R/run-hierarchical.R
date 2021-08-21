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
library(ape)
library(dplyr)

options(future.globals.maxSize = 1024 * 1024^2)
sessionInfo()

set.seed(42)
args <- commandArgs(trailingOnly = TRUE)
scale_method <- args[1]
cluster_method <- args[2]

# %% read data
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
colnames(cluster_df) <- "clusters"
read_df <- read_df[, rownames(cluster_df)]
if (!all(rownames(cluster_df) == colnames(read_df))) {
    stop("AssertionError")
}

regions <- jsonlite::read_json(
    paste0(WORKDIR, "results/cluster/", scale_method, "-", cluster_method, "/regions.json"),
    simplifyVector = TRUE
    )$regions

# %% build mean_df
mean_df <- data.frame(row.names = rownames(read_df))
for (cluster in unique(cluster_df[, 1])) {
    subset_df <- filter(cluster_df, clusters == cluster)
    subset_df <- read_df[, rownames(subset_df)]
    mean_df[, cluster] <- rowMeans(subset_df)
}

# %% draw
regions_label <- list(
    cortex = 1,
    thalamus = 2,
    hypothalamus = 3,
    olfactory = 4,
    hippocampus = 5
)
tip_color <- c()
for (i in colnames(mean_df)) {
    flag <- 1
    for (j in names(regions)) {
        for (k in regions[[j]]) {
            if (i == k) {
                tip_color <- c(tip_color, regions_label[[j]])
                flag = 0
            }
        }
    }
    if (flag) {
        tip_color <- c(tip_color, 6)
    }
}
names(tip_color) <- colnames(mean_df)
colors = c("violet", "peru", "red", "blue", "green", "black")
hc <- hclust(dist(t(mean_df)))
jpeg(
    paste0(
        WORKDIR, "results/cluster/", scale_method, "-", cluster_method, "/hierarchical.jpg"
    ),
    width = 1500, height = 1500
)
par(cex = 1.5)
plot(as.phylo(hc), type = "fan", tip.color = colors[tip_color], underscore = TRUE)
par(cex = 2)
legend(x = "topright", legend = names(regions_label), text.col = colors)
dev.off()
