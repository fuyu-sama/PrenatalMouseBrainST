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
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from umap import UMAP

session_info.show()
plt.rcParams.update({"font.size": 24})
session_info.show()
ncs = 12

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
regions_label = dict(cortex=0, hippocampus=1, thalamus=2, hypothalamus=3)

# %% read data
count_path = Path.joinpath(
    Path.home(), "workspace/mouse-brain-full/logcpm/scale_df/full-logcpm-inter.csv"
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    Path.home(), f"workspace/mouse-brain-full/SCT/SC3/full-{ncs}-SC3.csv"
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
count_full_df = count_full_df.reindex(cluster_df.index)

regions = dict(
    cortex=(
        "E135A_1",
        "E135B_11",
        "E155A_3",
        "E155B_9",
        "E165A_1",
        "E165B_10",
        "E175A1_4",
        "E175A2_2",
        "E175B_6",
        "P0A1_6",
        "P0A2_1",
        "P0B_1",
        "P0B_3",
    ),  # 皮层
    hippocampus=(
        "E135A_9",
        "E135B_2",
        "E155A_2",
        "E155B_8",
        "E165A_2",
        "E165B_9",
        "E175A1_6",
        "E175A2_5",
        "E175A2_6",
        "E175B_9",
        "P0A1_9",
        "P0A2_6",
        "P0B_4",
    ),  # 海马
    thalamus=(
        "E135A_2",
        "E135B_10",
        "E155A_12",
        "E155B_10",
        "E165A_3",
        "E165B_3",
        "E175A1_1",
        "E175A1_11",
        "E175A2_4",
        "E175A2_11",
        "E175B_5",
        "P0A1_8",
        "P0A2_5",
        "P0B_6",
        "P0B_7",
    ),  # 丘脑
    hypothalamus=(
        "E135A_6",
        "E135A_7",
        "E135B_3",
        "E135B_9",
        "E155A_11",
        "E155B_2",
        "E165A_4",
        "E165A_12",
        "E165B_1",
        "E175A1_9",
        "E175A2_10",
        "E175B_4",
        "P0A1_2",
        "P0A2_3",
        "P0B_12",
    ),  # 下丘脑
)

in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"sc3_{ncs}_clusters"])) if i not in in_regions
]

draw_cluster = cluster_df[f"sc3_{ncs}_clusters"].copy()
for c in regions:
    draw_cluster.replace(regions[c], c, inplace=True)
flag = len(regions)
for c in others:
    draw_cluster.replace(c, flag, inplace=True)
for c in regions:
    draw_cluster.replace(c, regions_label[c], inplace=True)

# %% embed
model_embed = UMAP()
result_embed = model_embed.fit_transform(count_full_df)

# %% draw
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xticks([])
ax.set_yticks([])
ax.scatter(
    result_embed[:, 0],
    result_embed[:, 1],
    c=draw_cluster,
    cmap="rainbow",
    s=2,
)
fig.savefig("test.jpg")
