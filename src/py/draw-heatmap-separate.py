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
from math import ceil
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 20})

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
    idx = sys.argv[1]
    scale_method = sys.argv[2]
except IndexError:
    idx = "E175A1"
    scale_method = "logcpm-hotspot-6-logcpm-0.99-Ai-0_500"

# %% read counts
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/logcpm/{idx}-logcpm.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

with open(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster/{scale_method}/regions.json",
        )) as f:
    regions = json.load(f)
with open(Path.joinpath(
        WORKDIR,
        f"results/cluster/colors.json",
)) as f:
    colors = json.load(f)
chosen_regions = {
    "cortex": "#CC163A",
    "hippocampus": "#4DFF4D",
    "thalamus": "#6495ED",
    "hypothalamus-mid": "#EE82EE"
}

# %% draw
# draw_spots
draw_spots = pd.DataFrame()
ncs = len(regions[idx])
cluster_df = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/{scale_method}/{idx}-{ncs}/{idx}-spots.csv",
    ),
    index_col=0,
    header=0,
)[f"type_1"]
cluster_df.name = "region_tip"
chosen_clusters = {}
for i in chosen_regions:
    for j in regions[idx]:
        if regions[idx][j] == "hypothalamus-top":
            # chosen_clusters[j] = "thalamus"
            pass
        else:
            # if regions[idx][j].startswith(i):
            if regions[idx][j] == i:
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
gene_path = Path.joinpath(
    WORKDIR,
    f"results/gene-cluster/{scale_method}/{idx}-{ncs}/{idx}-genes.csv",
)
gene_df = pd.read_csv(gene_path, index_col=0, header=0)
draw_genes = pd.DataFrame()
for i in chosen_clusters:
    temp_df = gene_df[gene_df["0"] == int(i)]
    temp_df = temp_df.replace(int(i), chosen_clusters[i])
    draw_genes = pd.concat([draw_genes, temp_df])
for j, i in enumerate(chosen_regions):
    draw_genes = draw_genes.replace(i, j)
draw_genes = draw_genes.sort_values(by="0")

# x axis genes
wanted_genes = [
    [
        "Foxg1",
        "Neurod6",
        "Bhlhe22",
        "Id2",
        "Zeb2",
        "Zbtb18",
        "Satb2",
        "Mef2c",
    ],
    [
        "Zbtb20",
        "Sema3c",
        "Eomes",
        "Fbln2",
        "Insm1",
        "Sstr2",
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
        "Pmch",
        "Slc16a13",
        "Dcn",
    ],
]
draw_genes_list = []
wanted_genes = [[j for j in i if j in draw_genes.index] for i in wanted_genes]
for i, gene_set in enumerate(wanted_genes):
    cluster_i_gene = list(draw_genes[draw_genes["0"] == i].index)
    cluster_i_gene = [j for j in cluster_i_gene if j not in gene_set]
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
rect_heatmap = [
    left + (spacing + cluster_width) * 2,
    bottom,
    width,
    height,
]

fig = plt.figure(figsize=(30, 10))
ax_heatmap = fig.add_axes(rect_heatmap)
ax_cluster = fig.add_axes(rect_cluster)

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
ax_cluster.set_yticklabels(
    ["neocortex", "hippocampus", "thalamus", "hypothalamus"])

[axis.set_visible(False) for axis in ax_heatmap.spines.values()]
[axis.set_visible(False) for axis in ax_cluster.spines.values()]

fig.savefig(
    Path.joinpath(WORKDIR, f"results/heatmap.{idx}.jpg"),
    bbox_inches="tight",
)
plt.close(fig)
