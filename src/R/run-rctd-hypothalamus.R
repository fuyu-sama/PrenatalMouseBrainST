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
library(RCTD)
library(Seurat)
library(dplyr)

sessionInfo()

set.seed(42)

sc_list <- list(
    E155A = "GSM3890934_Hypothalamus_traject_E15_srt_annotated_wo_blood.rds",
    E155B = "GSM3890934_Hypothalamus_traject_E15_srt_annotated_wo_blood.rds",
    E175A1 = "GSM3890935_Hypothalamus_traject_E17_srt_annotated_wo_blood.rds",
    E175A2 = "GSM3890935_Hypothalamus_traject_E17_srt_annotated_wo_blood.rds",
    E175B = "GSM3890935_Hypothalamus_traject_E17_srt_annotated_wo_blood.rds",
    P0A1 = "GSM3890936_Hypothalamus_traject_P0_srt_annotated_wo_blood.rds",
    P0A2 = "GSM3890936_Hypothalamus_traject_P0_srt_annotated_wo_blood.rds"
)

args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
scale_method <- args[2]
cluster_method <- args[3]
if (is.na(idx) | is.na(scale_method) | is.na(cluster_method)) {
    idx <- "E175B"
    scale_method <- "combat"
    cluster_method <- "sc3"
}

# %% read sc data and build reference
if (!(idx %in% names(sc_list))) quit(save = "no", status = 404)
seurat_obj <- readRDS(
    paste0(
        Sys.getenv("HOME"),
        "/Data/scRNAseq/2020_Nature_MolecularDesignMouseHypothalamus/COUNT/",
        sc_list[[idx]]
    )
)
sc_df <- as.data.frame(GetAssayData(seurat_obj, slot = "counts"))

cell_types <- as.character(Idents(seurat_obj))
for (i in 1:length(cell_types)) {
    cell_types[i] <- gsub("/", "or", cell_types[i])
}
names(cell_types) <- colnames(sc_df)
cell_types <- as.factor(cell_types)

removed_cells <- c()
for (i in unique(cell_types)) {
    if (length(cell_types[cell_types == i]) < 25) {
        removed_cells <- c(removed_cells, names(cell_types[cell_types == i]))
    }
}
cell_types <- cell_types[!(names(cell_types) %in% removed_cells)]
sc_df <- sc_df[, names(cell_types)]

n_umi_sc <- colSums(sc_df)
names(n_umi_sc) <- colnames(sc_df)

reference <- Reference(sc_df, cell_types, n_umi_sc)

# %% read st data
st_df <- read.csv(
    paste0(WORKDIR, "Data/scale_df/raw/", idx, "-raw.csv"),
    check.names = FALSE, row.names = 1
)

coor_df <- read.csv(
    paste0(WORKDIR, "Data/coor_df/", idx, "-coor.csv"),
    check.names = FALSE, row.names = 1
)

cluster_df <- read.csv(
    paste0(
        WORKDIR,
        "results/cluster/", scale_method, "-", cluster_method,
        "/pattern/full-", cluster_method, ".csv"
        ),
    check.names = F, row.names = 1
)
cluster_df <- as.data.frame(cluster_df[colnames(st_df), ])
rownames(cluster_df) <- colnames(st_df)
colnames(cluster_df) <- "clusters"

regions <- jsonlite::read_json(
    paste0(WORKDIR, "results/cluster/", scale_method, "-", cluster_method, "/regions.json"),
    simplifyVector = TRUE
    )$regions

cluster_df <- filter(
    cluster_df,
    clusters %in% regions[["hypothalamus"]]
)
st_df <- st_df[, rownames(cluster_df)]
coor_df <- coor_df[rownames(cluster_df), ]

n_umi_st <- colSums(st_df)
puck <- SpatialRNA(coor_df, st_df, n_umi_st)

# %% run RCTD
rctd_obj <- create.RCTD(puck, reference, max_cores = 10)
rctd_obj <- run.RCTD(rctd_obj, doublet_mode = "doublet")
for (i in names(rctd_obj@results)) {
    if (i == "score_mat") next
    write.csv(
        rctd_obj@results[[i]],
        paste0(
            WORKDIR, "results/RCTD/hypothalamus/", idx, "/results/", i, ".csv"
        )
    )
}
