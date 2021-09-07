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

# %% envoronment config
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

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
    scale_method = "combat"
    cluster_method = "sc3"

# %% read counts
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
count_full_df = count_full_df.reindex(cluster_df.index)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]
regions_label = {j: i for i, j in enumerate(regions)}
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"{cluster_method}_clusters"]))
    if i not in in_regions
]

# %% draw 1vsa
draw_cluster = cluster_df[f"{cluster_method}_clusters"].copy()
for c in regions:
    draw_cluster.replace(regions[c], c, inplace=True)
flag = len(regions_label)
for c in others:
    draw_cluster.replace(c, flag, inplace=True)
for c, i in zip(regions, range(len(regions))):
    draw_cluster.replace(c, regions_label[c], inplace=True)

for region in regions:
    genes = []  # up, avg_log2FC > 0
    for timepoint in timepoints:
        de_path = Path.joinpath(
            WORKDIR,
            f"results/DE/{scale_method}-{cluster_method}/timepoint-specific/DE-{region}-{timepoint}.csv",
        )
        try:
            de_df = pd.read_csv(de_path, index_col=0, header=0)
        except FileNotFoundError:
            continue
        de_df = de_df[(de_df["avg_log2FC"] > 0) & (de_df["p_val_adj"] <= 0.01)]
        de_df = de_df.sort_values(by="avg_log2FC", ascending=False)
        de_df.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/DE/{scale_method}-{cluster_method}/timepoint-specific/UP-{region}-{timepoint}.csv"
            ))
        for j in de_df.index:
            if j not in genes:
                genes.append(j)

    # build draw_df
    spots = draw_cluster[draw_cluster == regions_label[region]].sort_index()
    draw_df = count_full_df.reindex(spots.index, axis="index")
    draw_df = draw_df.reindex(genes, axis="columns")

    # pcolor
    flag = -1
    length = 0
    lastlength = 0
    lastlabel = ""
    timepoint_ticks = []
    timepoint_ticklabels = []
    pcolor_list = [i.split("_")[0] for i in spots.index]
    for i in range(len(pcolor_list)):
        if pcolor_list[i] != lastlabel:
            timepoint_ticks.append(int(length / 2) + lastlength)
            timepoint_ticklabels.append(lastlabel)
            lastlength += length
            length = 0
            flag += 1
        lastlabel = pcolor_list[i]
        length += 1
        pcolor_list[i] = flag
    timepoint_ticks.append(int(length / 2) + lastlength)
    timepoint_ticklabels.append(lastlabel)

    # z-score
    draw_df = (draw_df - draw_df.mean()) / draw_df.std()

    # draw
    left, bottom = 0.1, 0.1
    width, height = 0.65, 0.65
    spacing, cluster_width = 0.02, 0.01
    rect_heatmap = [left + spacing, bottom + height + spacing, width, height]
    rect_cluster = [left, bottom + height + spacing, cluster_width, height]
    fig = plt.figure(figsize=(15, 10))
    ax_heatmap = fig.add_axes(rect_heatmap)
    ax_cluster = fig.add_axes(
        rect_cluster,
        sharey=ax_heatmap,
    )
    hm = ax_heatmap.imshow(draw_df, cmap="bwr", aspect="auto", vmin=-3, vmax=3)
    plt.setp(ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
    ax_heatmap.set_xlabel(f"DEG in {region}: {draw_df.shape}")
    ax_heatmap.set_xticks([])
    ax_heatmap.xaxis.set_label_position("top")
    ax_cluster.pcolor(
        np.array(pcolor_list).reshape([len(pcolor_list), 1]),
        cmap=ListedColormap(colors[:len(regions) + len(others)]),
    )
    ax_cluster.set_xticks([])
    ax_cluster.set_yticks(timepoint_ticks)
    ax_cluster.set_yticklabels(timepoint_ticklabels)
    ax_cluster.set_xlabel("Regions")
    ax_cluster.xaxis.set_label_position("top")
    cb = fig.colorbar(hm, ax=ax_heatmap)

    [tk.set_visible(False) for tk in ax_heatmap.get_yticklabels()]
    [axis.set_visible(False) for axis in ax_heatmap.spines.values()]
    [axis.set_visible(False) for axis in ax_cluster.spines.values()]

    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/DE/{scale_method}-{cluster_method}/timepoint-specific/heatmap-{region}.jpg",
        ),
        bbox_inches="tight",
    )
    plt.close(fig)
