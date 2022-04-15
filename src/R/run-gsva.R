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
library(GSVA)
library(msigdb)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)
set.seed(42)
sessionInfo()

idx <- args[1]
scale_method <- args[2]

if (is.na(idx) | is.na(scale_method)) {
    idx <- "E165A"
    scale_method <- "logcpm"
}

# %% read data
count_df <- read.csv(
    paste0(
        WORKDIR,
        "Data/scale_df/", scale_method, "/", idx, "-", scale_method, ".csv"),
    check.names = F, row.names = 1
)

# %% get Msigdb
load_from_file <- FALSE
if (load_from_file) {
    load(paste0(WORKDIR, "Data/msigdb.mmu.rda"))
} else {
    msigdb.mmu <- getMsigdb(org = "mm", id = "SYM", version = "7.4")
    save(msigdb.mmu, file=paste0(WORKDIR, "Data/msigdb.mmu.rda"))
}

# %%
msigdb.mmu.ha <- msigdb.mmu[names(msigdb.mmu)[grepl('HALLMARK', names(msigdb.mmu))]]
myhallmark <- data.frame(set1 = msigdb.mmu.ha[[names(msigdb.mmu.ha)[1]]]@geneIds) %>%
    t(.) %>% data.frame(., check.names = FALSE)
for (i in names(msigdb.mmu.ha)[-1]){
    tmp <- data.frame(msigdb.mmu.ha[[i]]@geneIds) %>%
        t(.) %>% data.frame(., check.names = FALSE)
    myhallmark <- rbind.fill(myhallmark, tmp)
}
rownames(myhallmark) <- names(msigdb.mmu.ha)
write.table(
    myhallmark, paste0(WORKDIR, "Data/mouse_hallmark.gmt"),
    sep = '\t', quote = FALSE, row.names = TRUE, col.names = FALSE, na = ""
)

# %% perform GSVA
hallmark_gmt <- getGmt(paste0(WORKDIR, "Data/mouse_hallmark.gmt"))
gsva.es <- gsva(
    as.matrix(count_df),
    hallmark_gmt,
    parallel.sz = 10,
    verbose = FALSE
)
gsva.es <- rbind(id = colnames(gsva.es), gsva.es)
write.table(
    gsva.es, paste0(WORKDIR, "results/gsva/", idx, "-hallmark.csv"),
    quote = FALSE, col.names = FALSE
)
