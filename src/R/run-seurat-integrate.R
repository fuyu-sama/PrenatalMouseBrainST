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

options(future.globals.maxSize = 1024 * 1024^2)
sessionInfo()

idx_full <- list(
    E135A = "V10M17-100-E135A",
    E135B = "V10M17-085-E135B",
    E155A = "V10M17-100-E155A",
    E155B = "V10M17-085-E155B",
    E165A = "V10M17-100-E165A",
    E165B = "V10M17-085-E165B",
    E175A1 = "V10M17-101-E175A1",
    E175A2 = "V10M17-101-E175A2",
    E175B = "V10M17-085-E175B",
    P0B = "V10M17-100-P0B",
    P0A1 = "V10M17-101-P0A1",
    P0A2 = "V10M17-101-P0A2"
)

# %% read data
read_df_list <- list()
select_genes <- c()
for (idx in names(idx_full)) {
    read_df <- read.csv(
        paste0(WORKDIR, "Data/scale_df/raw_count/", idx, "-raw.csv"),
        check.names = FALSE, row.names = 1
    )
    if (length(select_genes) == 0) {
        select_genes <- rownames(read_df)
    } else {
        select_genes <- intersect(rownames(read_df), select_genes)
    }
    read_df_list[[idx]] <- read_df
}

# %% build seurat object list
seurat_list <- list()
for (idx in names(read_df_list)) {
    seurat_obj <- CreateSeuratObject(read_df_list[[idx]][select_genes, ])
    seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    seurat_obj <- FindVariableFeatures(
        seurat_obj, selection.method = "vst", nfeatures = 12000, verbose = FALSE
    )
    seurat_list[[idx]] <- seurat_obj
}

# %% perform integration
anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    anchor.features = length(select_genes),
    verbose = FALSE
)
seurat_combined <- IntegrateData(anchorset = anchors, verbose = FALSE)

write.csv(
    as.data.frame(seurat_combined[["integrated"]]@data),
    paste0(
        WORKDIR,
        "Data/scale_df/seurat_integrate/full-seurat_integrate-inter.csv"
    )
)
