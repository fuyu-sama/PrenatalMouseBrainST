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
from pathlib import Path

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

# %% read counts
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/logcpm/full-logcpm-inter.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/SCT-SC3/pattern/full-SC3.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
count_full_df = count_full_df.reindex(cluster_df.index)

regions = dict(
    cortex=("E135A_1", "E135B_3", "E135B_9", "E155A_4", "E155B_5", "E155B_6",
            "E165A_1", "E165B_4", "E175A1_8", "E175A2_5", "E175A2_6",
            "E175B_6", "P0A1_5", "P0A2_1", "P0B_3"),
    thalamus=("E135A_2", "E135B_6", "E135B_7", "E155A_5", "E155A_7", "E155A_8",
              "E155B_7", "E165A_3", "E165A_5", "E165B_3", "E175A1_2",
              "E175A1_7", "E175A2_10", "E175A2_11", "E175A2_12", "E175B_2",
              "E175B_4", "E175B_5", "P0A1_6", "P0A1_7", "P0A1_12", "P0A2_6",
              "P0B_1", "P0B_7", "P0B_8"),
    hypothalamus=("E135A_3", "E135A_6", "E135B_5", "E135B_12", "E155A_6",
                  "E155B_1", "E165A_6", "E165B_1", "E175B_3", "E175A1_1",
                  "E175A2_1", "P0A1_1", "P0A2_4", "P0B_12"),
    olfactory=("E135A_5", "E135A_7", "E135B_4", "E155A_11", "E155B_2",
               "E155B_11", "E165A_10", "E165A_11", "E165B_7", "E165B_10",
               "E175A1_5", "E175A2_7", "E175A2_9", "E175B_1", "E175B_9",
               "P0A1_10", "P0A2_7", "P0B_5", "P0B_6"),
    hippocampus=("E155A_2", "E165A_2", "E165B_6", "E175A2_14", "E175A2_15",
                 "E175B_11", "P0A1_8", "P0A2_2", "P0B_1"),
)
regions_label = dict(cortex=0,
                     thalamus=1,
                     hypothalamus=2,
                     olfactory=3,
                     hippocampus=4)
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"sc3_clusters"])) if i not in in_regions
]

# %% draw 1vsa
genes = []  # up, avg_log2FC > 0
# in regions
for region in regions:
    de_path = Path.joinpath(
        WORKDIR,
        f"results/DE/region-specific/DE-{region}.csv",
    )
    de_df = pd.read_csv(de_path, index_col=0, header=0)
    de_df = de_df[(de_df["avg_log2FC"] > 0) & (de_df["p_val_adj"] <= 0.01)]
    de_df = de_df.sort_values(by="avg_log2FC", ascending=False)
    de_df.to_csv(
        Path.joinpath(WORKDIR, f"results/DE/region-specific/UP-{region}.csv"))
    for j in de_df.index:
        if j not in genes:
            genes.append(j)

# x_ticks genes
gene_ticks = []
gene_tick_labels = []
wanted = ["Neurog2", "Sox11", "Tmem132b", "Otp", "Sox2", "Zbtb20", "Satb2"]
for i in range(len(genes)):
    if genes[i] in wanted:
        gene_ticks.append(i)
        gene_tick_labels.append(genes[i])

# cluster pcolor
region_ticks = []
length = [0]
draw_cluster = cluster_df[f"sc3_clusters"].copy()
for c in regions:
    draw_cluster.replace(regions[c], c, inplace=True)
flag = len(regions)
for c in others:
    draw_cluster.replace(c, flag, inplace=True)
    # flag += 1
for c, i in zip(regions, range(len(regions))):
    c_len = draw_cluster[draw_cluster == c].shape[0]
    region_ticks.append(c_len / 2 + sum(length[:i + 1]))
    length.append(c_len)
    draw_cluster.replace(c, regions_label[c], inplace=True)
draw_cluster.sort_values(inplace=True)

# build draw_df
draw_df = count_full_df.reindex(draw_cluster.index, axis="index")
draw_df = draw_df.reindex(genes, axis="columns")

# z-score
draw_df = (draw_df - draw_df.mean()) / draw_df.std()

# draw
left, bottom = 0.1, 0.1
width, height = 0.65, 0.65
spacing, cluster_width = 0.02, 0.01
rect_heatmap = [left + spacing, bottom + height + spacing, width, height]
rect_cluster = [left, bottom + height + spacing, cluster_width, height]
fig = plt.figure(figsize=(25, 10))
ax_heatmap = fig.add_axes(rect_heatmap)
ax_cluster = fig.add_axes(
    rect_cluster,
    sharey=ax_heatmap,
)
hm = ax_heatmap.imshow(draw_df, cmap="bwr", aspect="auto", vmin=-3, vmax=3)
plt.setp(ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
ax_heatmap.set_xlabel("Differential Expression Genes")
ax_heatmap.set_xticks(gene_ticks)
ax_heatmap.set_xticklabels(gene_tick_labels)
ax_heatmap.xaxis.set_label_position("top")
ax_cluster.pcolor(
    draw_cluster.to_numpy().reshape([len(draw_cluster), 1]),
    cmap=ListedColormap(colors[:len(regions) + len(others)]),
)
ax_cluster.set_xticks([])
ax_cluster.set_yticks(region_ticks)
ax_cluster.set_yticklabels(regions.keys())
ax_cluster.set_xlabel("Regions")
ax_cluster.xaxis.set_label_position("top")
cb = fig.colorbar(hm, ax=ax_heatmap)

[tk.set_visible(False) for tk in ax_heatmap.get_yticklabels()]
[axis.set_visible(False) for axis in ax_heatmap.spines.values()]
[axis.set_visible(False) for axis in ax_cluster.spines.values()]

fig.savefig(
    Path.joinpath(WORKDIR, "results/DE/region-specific/heatmap.jpg"),
    bbox_inches="tight",
)
plt.close(fig)
