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
plt.rcParams.update({"font.size": 12})

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
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
try:
    scale_method = sys.argv[1]
except IndexError:
    scale_method = "combat-Ai-500union"

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

with open(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/color_swift.json",
        )) as f:
    color_swift = json.load(f)
color_swift.pop("E135A", None)
chosen_colors = ["#CC163A", "#FFD700", "#EE82EE"]

# %% draw
# draw_spots
draw_spots = pd.DataFrame()
for idx in color_swift:
    cluster_df = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/pattern/{idx}-sc3.csv",
        ),
        index_col=0,
        header=0,
    )["sc3_12_clusters"]
    cluster_df.name = "region_tip"
    chosen_clusters = []
    for i in chosen_colors:
        for j in color_swift[idx]:
            if color_swift[idx][j] == i:
                chosen_clusters.append(j)
    for i, j in zip(chosen_clusters, chosen_colors):
        i = int(i)
        temp_df = cluster_df[cluster_df == i]
        temp_df = temp_df.replace(i, j).to_frame()
        temp_df["idx"] = [idx] * len(temp_df)
        draw_spots = pd.concat([draw_spots, temp_df])

for j, i in enumerate(chosen_colors):
    draw_spots = draw_spots.replace(i, j)
for j, i in enumerate(color_swift):
    draw_spots = draw_spots.replace(i, j)
draw_spots.sort_values(by=["region_tip", "idx"], inplace=True)

region_ticks = []
j = 0
for i in set(draw_spots["region_tip"]):
    ticks_len = len(draw_spots[draw_spots["region_tip"] == i])
    region_ticks.append(j + ticks_len / 2)
    j += ticks_len

# draw_genes
with open(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster/logcpm-hotspot-6-Ai-500union/cluster.json",
        )) as f:
    gene_cluster_json = json.load(f)
init_dict = {"cortex": 0, "thalamus": 0, "hypothalamus": 0}
gene_dict = {i: init_dict.copy() for i in count_full_df.columns}
for idx in color_swift:
    gene_path = Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/logcpm-hotspot-6-Ai-500union",
        f"{idx}-{gene_cluster_json[idx]['ncs']}/{idx}-genes.csv",
    )
    gene_df = pd.read_csv(gene_path, index_col=0, header=0)
    gene_cluster_region = gene_cluster_json[idx]["cluster"]
    for gene in gene_df.index:
        for region in gene_cluster_region:
            if gene_df.loc[gene, "0"] == gene_cluster_region[region]:
                gene_dict[gene][region] += 1

draw_genes = pd.Series(name="region")
indel = {"cortex": 1, "thalamus": 2, "hypothalamus": 3}
for gene in gene_dict:
    for region in gene_dict[gene]:
        if gene_dict[gene][region] >= len(color_swift) / 2:
            draw_genes[gene] = indel[region]
draw_genes.sort_values(inplace=True)

# build draw_df
draw_df = count_full_df.reindex(
    index=draw_spots.index,
    columns=draw_genes.index,
)

# z-score
draw_df = (draw_df - draw_df.mean()) / draw_df.std()

# draw
left, bottom = 0.1, 0.1
width, height = 0.65, 0.65
spacing, cluster_width = 0.01, 0.01
tip_width, tip_height = 0.2, 0.02
rect_cluster = [
    left,
    bottom,
    cluster_width,
    height,
]
rect_idx = [
    left + spacing + cluster_width,
    bottom,
    cluster_width,
    height,
]
rect_heatmap = [
    left + (spacing + cluster_width) * 2,
    bottom,
    width,
    height,
]
rect_idx_tip = [
    left + spacing * 2 + cluster_width,
    bottom + height + spacing,
    tip_width,
    tip_height,
]

fig = plt.figure(figsize=(20, 10))
ax_heatmap = fig.add_axes(rect_heatmap)
ax_cluster = fig.add_axes(rect_cluster)
ax_idx = fig.add_axes(rect_idx)
ax_idx_tip = fig.add_axes(rect_idx_tip)

hm = ax_heatmap.imshow(draw_df, cmap="bwr", aspect="auto", vmin=-2, vmax=2)
plt.setp(ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
ax_heatmap.xaxis.set_label_position("top")
ax_heatmap.set_xticks([])
ax_heatmap.set_yticks([])
cb = fig.colorbar(hm, ax=ax_heatmap)

ax_cluster.pcolor(
    draw_spots["region_tip"].to_numpy().reshape([len(draw_spots), 1]),
    cmap=ListedColormap(chosen_colors),
)
ax_cluster.invert_yaxis()
ax_cluster.set_xticks([])
ax_cluster.set_yticks(region_ticks)
# ax_cluster.set_yticklabels(["cortex", "ventricle", "thalamus", "hypothalamus"])
ax_cluster.set_yticklabels(["cortex", "thalamus", "hypothalamus"])

ax_idx.pcolor(
    draw_spots["idx"].to_numpy().reshape([len(draw_spots), 1]),
    cmap="Paired",
)
ax_idx.invert_yaxis()
ax_idx.set_xticks([])
ax_idx.set_yticks([])

idx_tip = np.array([i for i in range(len(color_swift))])
ax_idx_tip.pcolor(
    idx_tip.reshape([1, len(idx_tip)]),
    cmap="Paired",
)
ax_idx_tip.set_xticks([i + 0.5 for i in range(len(color_swift))])
ax_idx_tip.set_xticklabels(color_swift.keys())
ax_idx_tip.set_yticks([])
ax_idx_tip.xaxis.tick_top()
plt.setp(ax_idx_tip.get_xticklabels(), rotation=45)

[axis.set_visible(False) for axis in ax_heatmap.spines.values()]
[axis.set_visible(False) for axis in ax_cluster.spines.values()]
[axis.set_visible(False) for axis in ax_idx.spines.values()]
[axis.set_visible(False) for axis in ax_idx_tip.spines.values()]

fig.savefig(
    Path.joinpath(WORKDIR, f"results/heatmap.jpg"),
    bbox_inches="tight",
)
plt.close(fig)
