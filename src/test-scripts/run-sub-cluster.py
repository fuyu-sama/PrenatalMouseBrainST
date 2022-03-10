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

# %%
import json
import os
import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image
from sklearn.cluster import KMeans

from SpaGene.spagene import Spagene

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
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
    "E175B": "V10M17-085-E175B",
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
    cluster_method = sys.argv[2]
    idx = sys.argv[3]
    region = sys.argv[4]

except IndexError:
    scale_method = "combat-zjq"
    cluster_method = "sc3"
    idx = "E165A"
    region = "hypothalamus"

# %%
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
)
count_df_all = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

coor_path = Path.joinpath(
    WORKDIR,
    f"Data/coor_df/{idx}-coor.csv",
)
coor_df_all = pd.read_csv(
    coor_path,
    index_col=0,
    header=0,
)

de_path = Path.joinpath(
    WORKDIR,
    f"results/DE/{scale_method}-{cluster_method}/region-specific/DE-{region}.csv",
)
de_df = pd.read_csv(de_path, index_col=0, header=0)
de_df = de_df[(de_df["avg_log2FC"] > 0) & (de_df["p_val_adj"] <= 0.01)]
de_df = de_df.sort_values(by="avg_log2FC", ascending=False)

he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = Image.open(he_path)

cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

in_regions = [
    i for i in count_df_all.index
    if cluster_df.loc[i, "sc3_clusters"] in regions[region]
]

count_df = count_df_all.reindex(index=in_regions, columns=de_df.index)
coor_df = coor_df_all.reindex(index=count_df.index)

# %%

# %%
cluster_result = pd.DataFrame(index=coor_df.index)
count_df_1 = count_df.reindex(columns=spagenes.index)
fig, axes = plt.subplots(1, 5, figsize=(50, 10))
[ax.axis("off") for ax in axes.flatten()]
[ax.imshow(he_image) for ax in axes.flatten()]
for n, ax in enumerate(axes.flatten()):
    model_cluster = KMeans(n_clusters=n + 2)
    cluster_result[n + 2] = model_cluster.fit_predict(count_df_1)
    ax.set_title(f"ncs = {n + 2}")
    for i in range(n + 2):
        draw_df = coor_df.reindex(index=[
            j for j in cluster_result[n + 2].index
            if i == cluster_result.loc[j, n + 2]
        ], )
        ax.scatter(
            draw_df["X"],
            draw_df["Y"],
            c=colors[i],
            s=16,
            label=i,
        )
    ax.legend()

fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-cluster.jpg"
    ))
count_df_1.T.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-count.csv"
    ))
cluster_result.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-cluster.csv"
    ))

# %%
for ncs in range(2, n + 3):
    os.system(
        f"Rscript {Path.joinpath(WORKDIR, 'src/R/run-sub-de.R')} {scale_method} {cluster_method} {idx} {region} {ncs}"
    )

    de_df_1 = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-{ncs}-de.csv"
        ),
        index_col=0,
        header=0,
    )
    de_df_1 = de_df_1[(de_df_1["avg_log2FC"] > 0) & (de_df_1["p_val_adj"] <= 0.01)]
    de_df_1 = de_df_1.sort_values(
        by=["cluster", "avg_log2FC"],
        ascending=[True, False],
    )

    draw_cluster = cluster_result[ncs].sort_values(ascending=True)

    draw_genes = []
    [draw_genes.append(i) for i in de_df_1["gene"] if i not in draw_genes]

    draw_df = count_df_1.reindex(
        columns=draw_genes,
        index=draw_cluster.index,
    )
    draw_df = (draw_df - draw_df.mean()) / draw_df.std()

    left, bottom = 0.1, 0.1
    width, height = 0.8, 0.4
    spacing, cluster_width = 0.02, 0.01
    rect_he = [left, bottom + height + spacing, height, height]
    rect_heatmap = [
        left + height * 0.9 + spacing * 2, bottom + height + spacing, width,
        height
    ]
    rect_cluster = [
        left + height * 0.9 + spacing * 1, bottom + height + spacing,
        cluster_width, height
    ]
    fig = plt.figure(figsize=(30, 25))
    ax_he = fig.add_axes(rect_he)
    ax_heatmap = fig.add_axes(rect_heatmap)
    ax_cluster = fig.add_axes(
        rect_cluster,
        sharey=ax_heatmap,
    )

    ax_he.imshow(he_image)
    ax_he.scatter(
        coor_df["X"].reindex(index=draw_cluster.index),
        coor_df["Y"].reindex(index=draw_cluster.index),
        s=16,
        c=draw_cluster,
        cmap=ListedColormap(colors[:ncs]),
    )
    ax_he.set_title(f"n_clusters = {ncs}")
    ax_he.axis("off")

    hm = ax_heatmap.imshow(draw_df, cmap="bwr", aspect="auto", vmin=-3, vmax=3)
    plt.setp(ax_heatmap.get_xticklabels(), rotation=90, fontsize=10)
    ax_heatmap.set_xticks(range(draw_df.shape[1]))
    ax_heatmap.set_xticklabels(draw_df.columns)
    ax_heatmap.xaxis.set_label_position("top")
    ax_heatmap.set_title(
        f"de genes({count_df.shape[1]}) -> spatial genes({spagenes.shape[0]}) -> de genes({draw_df.shape[1]})"
    )
    ax_cluster.pcolor(
        draw_cluster.to_numpy().reshape([len(draw_cluster), 1]),
        cmap=ListedColormap(colors[:ncs]),
    )
    ax_cluster.set_xticks([])
    ax_cluster.set_yticks([])
    ax_cluster.xaxis.set_label_position("top")
    cb = fig.colorbar(hm, ax=ax_heatmap)

    [tk.set_visible(False) for tk in ax_heatmap.get_yticklabels()]
    [axis.set_visible(False) for axis in ax_heatmap.spines.values()]
    [axis.set_visible(False) for axis in ax_cluster.spines.values()]

    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/{idx}-{region}-{ncs}.jpg"
        ),
        bbox_inches="tight",
    )
