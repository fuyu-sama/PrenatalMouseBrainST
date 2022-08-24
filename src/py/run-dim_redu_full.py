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
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 16})

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

# %% read cluster result
region_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-sc3/regions.json")
with open(region_path) as f:
    regions = json.load(f)
with open(Path.joinpath(WORKDIR, f"results/cluster/colors.json")) as f:
    colors = json.load(f)

count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/combat/full-combat.csv",
)
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T
chosen_genes = []
with open(Path.joinpath(WORKDIR, f"Data/gene-lists/{scale_method}.csv")) as f:
    for line in f:
        line = line.strip(line)
        chosen_genes.append(line)

clusters = pd.Series()
for idx in regions:
    ncs = len(regions[idx])
    cluster_path = Path.joinpath(
        WORKDIR,
        f"results/cluster/{scale_method}-sc3/pattern/{idx}-sc3.csv",
    )
    cluster_series = pd.read_csv(
        cluster_path,
        index_col=0,
        header=0,
    )[f"sc3_{ncs}_clusters"].astype(str)
    for i in regions[idx]:
        if regions[idx][i] in colors:
            flag = colors[regions[idx][i]]
        else:
            flag = regions[idx][i]
        cluster_series.replace(i, flag, inplace=True)
    clusters = pd.concat([clusters, cluster_series])

count_df = count_df.reindex(index=clusters.index, columns=chosen_genes)
chosen_regions = [
    "cortex", "hippocampus", "hypothalamus", "thalamus"
]

# %% tSNE
result_tsne = TSNE().fit_transform(
    PCA(n_components=40).fit_transform(count_df))
result_df = pd.DataFrame(
    result_tsne,
    index=count_df.index,
    columns=["X", "Y"],
)
result_df.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-sc3/full-tsne.csv",
    ))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
ax1.set_xticks([])
ax1.set_yticks([])
chosen_spots = []
for i in chosen_regions:
    draw_df = clusters[clusters == colors[i]]
    draw_df = result_df.reindex(index=draw_df.index)
    ax1.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=i,
        c=i,
        s=8,
    )
    [chosen_spots.append(i) for i in draw_df.index]
unchosen_spots = [
    i for i in clusters.index if clusters[i] not in chosen_spots
]
draw_df = result_df.reindex(unchosen_spots)
ax1.scatter(
    draw_df["X"],
    draw_df["Y"],
    c=clusters.reindex(index=draw_df.index),
    s=8,
)
ax1.legend(
    markerscale=2,
    loc='upper left',
    bbox_to_anchor=(-0.4, 1),
)

ax2.set_xticks([])
ax2.set_yticks([])
for i, idx in enumerate(regions):
    draw_df = result_df.reindex(index=[j for j in result_df.index if idx in j])
    ax2.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=idx_tp[idx],
        c=colors[i],
        s=8,
    )
ax2.legend(
    markerscale=2,
    loc='upper right',
    bbox_to_anchor=(1.3, 1),
)
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-sc3/full-tsne.jpg",
    ),
    bbox_inches="tight",
)
plt.close(fig)

# %% UMAP
result_umap = UMAP().fit_transform(
    PCA(n_components=40).fit_transform(count_df))
result_df = pd.DataFrame(
    result_umap,
    index=count_df.index,
    columns=["X", "Y"],
)
result_df.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-sc3/full-umap.csv",
    ))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
ax1.set_xticks([])
ax1.set_yticks([])
chosen_spots = []
for i in chosen_regions:
    draw_df = clusters[clusters == colors[i]]
    draw_df = result_df.reindex(index=draw_df.index)
    ax1.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=i,
        c=i,
        s=8,
    )
    [chosen_spots.append(i) for i in draw_df.index]
unchosen_spots = [
    i for i in clusters.index if clusters[i] not in chosen_spots
]
draw_df = result_df.reindex(unchosen_spots)
ax1.scatter(
    draw_df["X"],
    draw_df["Y"],
    c=clusters.reindex(index=draw_df.index),
    s=8,
)
ax1.legend(
    markerscale=2,
    loc='upper left',
    bbox_to_anchor=(-0.4, 1),
)

ax2.set_xticks([])
ax2.set_yticks([])
for i, idx in enumerate(regions):
    draw_df = result_df.reindex(index=[j for j in result_df.index if idx in j])
    ax2.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=idx_tp[idx],
        c=colors[i],
        s=8,
    )
ax2.legend(
    markerscale=2,
    loc='upper right',
    bbox_to_anchor=(1.3, 1),
)
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-sc3/full-umap.jpg",
    ),
    bbox_inches="tight",
)
plt.close(fig)
