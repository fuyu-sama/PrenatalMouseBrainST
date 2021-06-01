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

session_info.show()
plt.rcParams.update({"font.size": 24})
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
colors = [
    "#FAEBD7",
    "#00FFFF",
    "#FFD700",
    "#0000FF",
    "#FF8C00",
    "#EE82EE",
    "#9ACD32",
    "#5F9EA0",
    "#7FFF00",
    "#7FFFD4",
    "#6495ED",
    "#008B8B",
    "#B8860B",
    "#C0C0C0",
    "#000080",
    "#D8BFD8",
    "#00CED1",
    "#9400D3",
    "#8E804B",
    "#0089A7",
    "#CB1B45",
    "#FFB6C1",
    "#00FF00",
    "#800000",
    "#376B6D",
    "#D8BFD8",
    "#F5F5F5",
    "#D2691E",
]

# %% read counts
count_path = Path.joinpath(
    Path.home(),
    f"workspace/mouse-brain-full/logcpm/scale_df/full-logcpm-inter.csv")
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    Path.home(), f"workspace/mouse-brain-full/SCT/SC3/full-{ncs}-SC3.csv")
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
regions_label = dict(cortex=0, hippocampus=1, thalamus=2, hypothalamus=3)
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"sc3_{ncs}_clusters"]))
    if i not in in_regions
]

# %% draw 1vsa
genes = []  # up, avg_log2FC > 0
# in regions
for region in regions:
    de_path = Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/logcpm/DE/region-specific/DE-{ncs}-{region}.csv",
    )
    de_df = pd.read_csv(de_path, index_col=0, header=0)
    de_df = de_df[(de_df["avg_log2FC"] > 0) & (de_df["p_val_adj"] <= 0.01)]
    de_df = de_df.sort_values(by="avg_log2FC", ascending=False)
    de_df.to_csv(
        Path.joinpath(
            Path.home(),
            f"workspace/mouse-brain-full/logcpm/DE/region-specific/UP-{ncs}-{region}.csv"
        ))
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
draw_cluster = cluster_df[f"sc3_{ncs}_clusters"].copy()
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
    Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/logcpm/DE/region-specific/heatmap-{ncs}-zscore.jpg",
    ),
    bbox_inches="tight",
)
plt.close(fig)
