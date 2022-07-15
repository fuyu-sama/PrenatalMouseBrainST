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
#               \  .-\__  '-'  ___/-. /
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
library(SPARK)

args <- commandArgs(trailingOnly = TRUE)
set.seed(16)

idx <- args[1]

# %% read path
count_path <- paste0(WORKDIR, "Data/scale_df/raw/", idx, "-raw.csv")
coor_path <- paste0(WORKDIR, "Data/coor_df/", idx, "-coor.csv")

count_df <- read.csv(count_path, check.names = F, row.names = 1)
coor_df <- read.csv(coor_path, check.names = F, row.names = 1)

# %% find spatial pattern genes
spark_obj <- CreateSPARKObject(
    counts = count_df,
    location = coor_df,
)
spark_obj@lib_size <- apply(spark_obj@counts, 2, sum)
spark_obj <- spark.vc(
    spark_obj,
    covariates = NULL,
    lib_size = spark_obj@lib_size,
    num_core = 20,
    verbose = FALSE
)
spark_obj <- spark.test(
    spark_obj,
    check_positive = TRUE,
    verbose = FALSE
)
spark_result <- subset(
    spark_obj@res_mtest,
    adjusted_pvalue < 1,
    select = c("combined_pvalue", "adjusted_pvalue")
)
spark_result <- spark_result[order(spark_result, decreasing = FALSE), ]
spark_result <- na.omit(spark_result)
write.csv(
    x = spark_result,
    file = paste0(WORKDIR, "results/spark/", idx, "-spark.csv")
)
