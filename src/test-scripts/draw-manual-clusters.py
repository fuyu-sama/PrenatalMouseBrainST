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
from itertools import combinations
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from scipy.cluster import hierarchy as sch

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
tp_full = {
    "E135": ["E135A", "E135B"],
    "E155": ["E155A", "E155B"],
    "E165": ["E165A", "E165B"],
    "E175": ["E175A1", "E175A2"],
    # "E175": ["E175A1", "E175A2", "E175B"],
    "P0": ["P0A1", "P0A2"],
}
idx_tp = {
    "E135A": "E135",
    "E135B": "E135",
    "E155A": "E155",
    "E155B": "E155",
    "E165A": "E165",
    "E165B": "E165",
    "E175A1": "E175",
    "E175A2": "E175",
    # "E175B": "E175",
    # "P0B": "P0",
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
    scale_method = "logcpm-hotspot-6"

# %% read data
with open(Path.joinpath(WORKDIR, "Data/clusters.json")) as f:
    gene_clusters = json.load(f)

for idx in idx_full:
    count_path = Path.joinpath(
        WORKDIR,
        f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
    )
    count_df = pd.read_csv(
        count_path,
        index_col=0,
        header=0,
    ).T

    coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)

    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    he_image = Image.open(he_path)

    all_genes = set([j for i in gene_clusters.values() for j in i])
    count_df = count_df.reindex(columns=all_genes)

    # cluster
    n_gene_clusters = 3
    n_spot_clusters = 6

    gene_distmat = sch.distance.pdist(count_df.T, metric="jaccard")
    spot_distmat = sch.distance.pdist(count_df, metric="jaccard")

    Z_gene = sch.linkage(gene_distmat, method="ward")
    gene_result = pd.Series(
        sch.fcluster(Z_gene, t=n_gene_clusters, criterion="maxclust"),
        index=count_df.columns,
    ).sort_values()

    Z_spot = sch.linkage(spot_distmat, method="ward")
    spot_result = pd.Series(
        sch.fcluster(Z_spot, t=n_spot_clusters, criterion="maxclust"),
        index=count_df.index,
    ).sort_values()

    fig, ax_heatmap = plt.subplots(figsize=(10, 10), dpi=50)

    ax_heatmap.imshow(
        count_df.reindex(columns=gene_result.index, index=spot_result.index),
        cmap="Reds",
        aspect="auto",
    )
    ax_heatmap.set_title(f"{idx}")
    ax_heatmap.set_xticks(range(len(gene_result)))
    ax_heatmap.set_xticklabels(list(gene_result.index), rotation=90)
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("Genes")
    ax_heatmap.set_ylabel("Spots")

    flag = 0
    for i in range(1, n_gene_clusters):
        flag += len(gene_result[gene_result == i])
        ax_heatmap.plot([flag - 0.5, flag - 0.5], [0, count_df.shape[0] - 5])

    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster-manual/{idx}-cluster.jpg",
        ),
        bbox_inches="tight",
    )
    gene_result.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster-manual/{idx}-cluster.csv",
        ))
    plt.close(fig)

    # draw distribution
    flag = 0
    for n_cb in range(len(gene_clusters)):
        for cb in combinations(gene_clusters.keys(), n_cb + 1):
            # cb: (cluster_1, cluster_2, ...)
            all_genes = []
            for cluster in cb:
                [all_genes.append(i) for i in gene_clusters[cluster]]
            all_genes = set(all_genes)
            draw_df = count_df.reindex(columns=all_genes).T.mean()

            fig, ax = plt.subplots(figsize=(10, 13), dpi=50)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.imshow(he_image)
            ax.axis("off")
            sc = ax.scatter(
                coor_df["X"],
                coor_df["Y"],
                c=draw_df,
                cmap="autumn_r",
                vmin=0,
                s=16,
                alpha=0.7,
            )
            cb = [f"\n({i})" for i in cb]
            title = f"{idx}{''.join(cb)}"
            title = f"{idx}\nall" if len(cb) == len(gene_clusters) else title
            ax.set_title(title)
            plt.colorbar(sc, orientation="horizontal")
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/gene-cluster-manual/{idx}-{flag}.jpg",
                ))
            flag += 1
            plt.close(fig)

    he_image.close()
