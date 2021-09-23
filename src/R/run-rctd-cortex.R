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
library(dplyr)

sessionInfo()

set.seed(42)

sc_list <- list(
    E135A = "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5",
    E135B = "GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5",
    E155A = "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5",
    E155B = "GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5",
    E165A = "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
    E165B = "GSM4635077_E16_filtered_gene_bc_matrices_h5.h5",
    E175A1 = "GSM5277844_E17_5_filtered_feature_bc_matrix.h5",
    E175A2 = "GSM5277844_E17_5_filtered_feature_bc_matrix.h5",
    E175B = "GSM5277844_E17_5_filtered_feature_bc_matrix.h5",
    P0A1 = "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
    P0A2 = "GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5"
)
timepoint_list <- list(
    E135A = "E13_5_",
    E135B = "E13_5_",
    E155A = "E15_5_",
    E155B = "E15_5_",
    E165A = "E16_5_",
    E165B = "E16_5_",
    E175A1 = "E17_5_",
    E175A2 = "E17_5_",
    E175B = "E17_5_",
    P0A1 = "P1_S1_",
    P0A2 = "P1_S1_"
)

args <- commandArgs(trailingOnly = TRUE)
idx <- args[1]
scale_method <- args[2]
cluster_method <- args[3]
if (is.na(idx) | is.na(scale_method) | is.na(cluster_method)) {
    idx <- "E165A"
    scale_method <- "combat"
    cluster_method <- "sc3"
}

# %% read sc data and build reference
sc_df <- as.data.frame(
    Seurat::Read10X_h5(
        paste0(
            Sys.getenv("HOME"),
            "/Data/scRNAseq/2021_Nature_MolecularLogicMouseBrain/COUNT/",
            sc_list[[idx]]
        )
    )
)
cell_names <- c()
for (i in strsplit(colnames(sc_df), "-")) {
    cell_names <- c(cell_names, paste0(timepoint_list[[idx]], i[1]))
}
colnames(sc_df) <- cell_names

meta_df <- read.delim(
    paste0(
        Sys.getenv("HOME"),
        "/Data/scRNAseq/2021_Nature_MolecularLogicMouseBrain/META/",
        "metaData_scDevSC.txt"
        ),
    row.names = 1, check.names = FALSE
)
meta_names <- c()
for (i in rownames(meta_df)) {
    if (grepl("-", i, fixed = TRUE)) {
        meta_names <- c(meta_names, strsplit(i, "-")[[1]][1])
    } else {
        meta_names <- c(meta_names, i)
    }
}
rownames(meta_df) <- meta_names

cell_names <- intersect(cell_names, rownames(meta_df))
meta_df <- meta_df[cell_names, ]
sc_df <- sc_df[, cell_names]

cell_types <- meta_df$New_cellType
names(cell_types) <- rownames(meta_df)

removed_cells <- c()
for (i in unique(cell_types)) {
    if (length(cell_types[cell_types == i]) < 25) {
        removed_cells <- c(removed_cells, names(cell_types[cell_types == i]))
    }
}
cell_types <- cell_types[!(names(cell_types) %in% removed_cells)]
sc_df <- sc_df[, names(cell_types)]
cell_types <- as.factor(cell_types)

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

st_df <- st_df[, rownames(coor_df)]

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
            WORKDIR, "results/RCTD/cortex/", idx, "/results/", i, ".csv"
        )
    )
}
