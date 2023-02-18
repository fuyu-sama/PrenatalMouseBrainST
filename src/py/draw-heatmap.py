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
from copy import deepcopy
from math import ceil
from pathlib import Path

import numpy as np
import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 26})

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
idx_tp = {
    "E135A": "E13.5",
    "E135B": "E13.5",
    "E155A": "E15.5",
    "E155B": "E15.5",
    "E165A": "E16.5",
    "E165B": "E16.5",
    "E175A1": "E17.5",
    "E175A2": "E17.5",
    "P0A1": "P0",
    "P0A2": "P0",
}

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
            f"results/cluster/{scale_method}-sc3/regions.json",
        )) as f:
    regions = json.load(f)
with open(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/colors.json",
        )) as f:
    colors = json.load(f)
chosen_regions = {
    "cortex": "#CC163A",
    "thalamus": "#6495ED",
    "hypothalamus": "#D2691E"
}

# %% draw
# draw_spots
draw_spots = pd.DataFrame()
for idx in regions:
    ncs = len(regions[idx])
    cluster_df = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/pattern/{idx}-sc3.csv",
        ),
        index_col=0,
        header=0,
    )[f"sc3_{ncs}_clusters"]
    cluster_df.name = "region_tip"
    chosen_clusters = {}
    for i in chosen_regions:
        for j in regions[idx]:
            if i in regions[idx][j]:
                chosen_clusters[j] = i
    for i in chosen_clusters:
        temp_df = cluster_df[cluster_df == int(i)]
        temp_df = temp_df.replace(int(i), chosen_clusters[i]).to_frame()
        temp_df["idx"] = [idx] * len(temp_df)
        draw_spots = pd.concat([draw_spots, temp_df])

for j, i in enumerate(chosen_regions):
    draw_spots = draw_spots.replace(i, j)
for j, i in enumerate(regions):
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
gene_dict = {i: deepcopy(init_dict) for i in count_full_df.columns}
for idx in regions:
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
    max_region = max(gene_dict[gene], key=gene_dict[gene].get)
    if gene_dict[gene][max_region] >= 1:
        draw_genes[gene] = indel[max_region]
draw_genes.sort_values(inplace=True)

# x axis genes
wanted_genes = [
    [
        "Mef2c",
        "Bcl11a",
        "Sox5",
        "Hivep2",
        "Satb2",
        "Fezf2",
        "Neurod2",
        "Neurod6",
        "Zbtb18",
        "Tbr1",
    ],
    [
        "Zfp423",
        "Tcf7l2",
        "Zic1",
        "Zic3",
        "Zic4",
        "Zic5",
        "Gbx2",
        "Id4",
        "Slc18a2",
        "Clybl",
        "Lhx9",
        "Foxp2",
        "Sox2",
        "Calb2",
    ],
    [
        "Dlx1",
        "Otp",
        "Asb4",
        "Trh",
        "Ddc",
        "Dlk1",
        "Nkx2-2",
        "Magel2",
        "Nap1l5",
        "Peg10",
    ],
]
draw_genes_list = []
for i, gene_set in enumerate(wanted_genes):
    cluster_i_gene = list(draw_genes[draw_genes == i + 1].index)
    [cluster_i_gene.remove(j) for j in gene_set if j in draw_genes.index]
    space_len = ceil(len(cluster_i_gene) / (len(gene_set) + 1))
    j = 0
    k = 0
    while j < len(cluster_i_gene):
        draw_genes_list.append(cluster_i_gene[j])
        if j % space_len == 0 and j != 0:
            draw_genes_list.append(gene_set[k])
            k += 1
        j += 1
    ll = k - len(gene_set)
    if ll != 0:
        for i in range(1, ll + 1):
            draw_genes_list.append(gene_set[-ll])

draw_genes = draw_genes.reindex(index=draw_genes_list)
xticks = []
xticklabels = []
wanted_genes_flatten = [j for i in wanted_genes for j in i]
for i, gene in enumerate(draw_genes.index):
    if gene in wanted_genes_flatten:
        xticks.append(i)
        xticklabels.append(gene)
draw_genes.to_csv(Path.joinpath(WORKDIR, f"results/heatmap-genes.csv"))

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
    bottom + height + spacing * 4,
    tip_width,
    tip_height,
]

fig = plt.figure(figsize=(40, 10))
ax_heatmap = fig.add_axes(rect_heatmap)
ax_cluster = fig.add_axes(rect_cluster)
ax_idx = fig.add_axes(rect_idx)
ax_idx_tip = fig.add_axes(rect_idx_tip)

hm = ax_heatmap.imshow(draw_df, cmap="bwr", aspect="auto", vmin=-2, vmax=2)
plt.setp(ax_heatmap.get_xticklabels(), rotation=90)
ax_heatmap.xaxis.set_label_position("top")
ax_heatmap.set_xticks(xticks)
ax_heatmap.set_xticklabels(xticklabels, fontstyle="italic")
ax_heatmap.set_yticks([])
cb = fig.colorbar(hm, ax=ax_heatmap)

ax_cluster.pcolor(
    draw_spots["region_tip"].to_numpy().reshape([len(draw_spots), 1]),
    cmap=ListedColormap(chosen_regions.values()),
)
ax_cluster.invert_yaxis()
ax_cluster.set_xticks([])
ax_cluster.set_yticks(region_ticks)
ax_cluster.set_yticklabels(["cortex\n+\nhippocampus", "thalamus", "hypothalamus"])

ax_idx.pcolor(
    draw_spots["idx"].to_numpy().reshape([len(draw_spots), 1]),
    cmap="Paired",
)
ax_idx.invert_yaxis()
ax_idx.set_xticks([])
ax_idx.set_yticks([])

idx_tip = np.array([i for i in range(len(regions))])
ax_idx_tip.pcolor(
    idx_tip.reshape([1, len(idx_tip)]),
    cmap="Paired",
)
ax_idx_tip.set_xticks([i + 0.5 for i in range(len(regions))])
ax_idx_tip.set_xticklabels([idx_tp[i] for i in regions.keys()])
ax_idx_tip.set_yticks([])
ax_idx_tip.xaxis.tick_top()

[axis.set_visible(False) for axis in ax_heatmap.spines.values()]
[axis.set_visible(False) for axis in ax_cluster.spines.values()]
[axis.set_visible(False) for axis in ax_idx.spines.values()]
[axis.set_visible(False) for axis in ax_idx_tip.spines.values()]

fig.savefig(
    Path.joinpath(WORKDIR, f"results/heatmap-{draw_df.shape[1]}_genes.svg"),
    bbox_inches="tight",
)
plt.close(fig)
