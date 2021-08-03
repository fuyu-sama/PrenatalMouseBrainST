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
library(ChIPseeker)

sessionInfo()
regions <- c("cortex", "thalamus", "hypothalamus", "olfactory", "hippocampus")

scale_method = "combat"

# %% read data
spompe <- GenomicFeatures::makeTxDbFromGFF(
    paste0(
        Sys.getenv("HOME"),
        "/refgenome/spaceranger/refdata-gex-mm10-2020-A/genes/genes.gtf")
)
peak_path <- paste0(WORKDIR, "Data/ENCODE/forebrain.bed")

de_list <- list()
for (region in regions) {
    de_path <- paste0(
        WORKDIR,
        "results/DE/", scale_method, "/region-specific/UP-", region, ".csv.out"
    )
    de_list[[region]] <- as.vector(t(read.csv(de_path, check.names = FALSE)))
}

# %% annotate peak
peak_anno <- annotatePeak(
    peak_path,
    TxDb = spompe, tssRegion = c(-1000, 1000), overlap = "all",
    addFlankGeneInfo = TRUE, flankDistance = 100000,
    genomicAnnotationPriority = c("Promoter"), verbose = FALSE
)

promoter_df <- as.data.frame(peak_anno@anno)
promoter_df <- promoter_df %>% filter(annotation == "Promoter")

# %% filter promoters
outs <- list()
for (region in regions) {
    de_promoter <- promoter_df %>% filter(geneId %in% de_list[[region]]) %>%
        select(seqnames, start, end, geneId)
    flag <- 0
    outs[[region]] <- data.frame()
    for (i in unique(de_promoter[, "geneId"])) {
        temp <- filter(de_promoter, geneId == i)
        chr <- as.character(temp[1, 1])
        start <- min(temp[, 2])
        end <- max(temp[, 3])
        outs[[region]] <- rbind(outs[[region]], c(chr, start, end, flag, i))
        flag <- flag + 1
    }
    colnames(outs[[region]]) <- c("chr", "start", "end", "peakId", "geneId")
    write.table(
        outs[[region]],
        paste0(
            WORKDIR,
            "results/DE/", scale_method,
            "/region-specific/motifInputs/", region, ".bed"
        ),
        quote = FALSE, sep = "\t", row.names = FALSE
    )
}
