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

import numpy as np
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
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
colors = [
    "#FAEBD7", "#00FFFF", "#FFD700", "#0000FF", "#FF8C00", "#EE82EE",
    "#9ACD32", "#5F9EA0", "#7FFF00", "#7FFFD4", "#6495ED", "#008B8B",
    "#B8860B", "#C0C0C0", "#000080", "#D8BFD8", "#00CED1", "#9400D3",
    "#8E804B", "#0089A7", "#CB1B45", "#FFB6C1", "#00FF00", "#800000",
    "#376B6D", "#D8BFD8", "#F5F5F5", "#D2691E"
]

timepoints = dict(
    E135=("E135A", "E135B"),
    E155=("E155A", "E155B"),
    E165=("E165A", "E165B"),
    E175=("E175A1", "E175A2", "E175B"),
    P0=("P0A1", "P0A2"),
)

try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
except IndexError:
    scale_method = "combat-qq-pattern"
    cluster_method = "sc3"

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

# %% build draw_df
region = "cortex"
genes = [
    "Arx", "Ascl1", "Dlx1", "Foxp4", "Gsx1", "Helt", "Six3", "Sp9"
]
draw_df = pd.DataFrame()
for i in timepoints:
    temp_df = count_df.reindex(index=[k for k in cluster_df.index if i in k])
    temp_df = temp_df.reindex(index=[
        k for k in temp_df.index
        if cluster_df.loc[k, "sc3_clusters"] in regions[region]
    ])
    for j in genes:
        if j not in count_df.columns:
            print(j)
            continue
        draw_df.loc[j, i] = np.mean(temp_df[j])

# %% draw
fig, ax = plt.subplots(figsize=(5, 10))
for x in range(draw_df.shape[1]):
    for y in range(draw_df.shape[0]):
        ax.scatter(x, y, s=draw_df.iloc[y, x] ** 2 * 100, c=colors[y + 1])
ax.set_xticks(range(draw_df.shape[1]))
ax.set_yticks(range(draw_df.shape[0]))
ax.set_xticklabels(draw_df.columns, rotation=45, ha="right")
ax.set_yticklabels(draw_df.index)
ax.set_title(region)
fig.savefig(Path.joinpath(WORKDIR, "results/test.jpg"), bbox_inches="tight")
