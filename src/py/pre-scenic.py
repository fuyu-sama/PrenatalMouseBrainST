#! venv/bin/python
# -*- encoding: utf-8 -*-

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
import json
import sys
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
    region = sys.argv[3]
except IndexError:
    scale_method = "combat-qq-logcpm"
    cluster_method = "sc3"
    region = "cortex"

# %% read count data
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
count_df = pd.read_csv(count_path, index_col=0, header=0).T
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

# %% do
mean_df = pd.DataFrame(index=count_df.columns)
for idx in idx_full:
    sub_df = cluster_df.reindex(index=[i for i in cluster_df.index if idx in i])
    for region in regions:
        in_region = []
        for i in sub_df.index:
            if sub_df.loc[i, "sc3_clusters"] in regions[region]:
                in_region.append(i)
        region_df = count_df.reindex(index=in_region)
        mean_df[f"{idx}_{region}"] = region_df.mean()

mean_df = mean_df.dropna(axis=1)
mean_df.T.to_csv(
    Path.joinpath(WORKDIR, f"results/SCENIC/{scale_method}-{cluster_method}/mean.csv")
)
