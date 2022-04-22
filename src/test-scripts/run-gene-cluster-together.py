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
from math import ceil
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from PIL import Image
from scipy.cluster import hierarchy as sch
from sklearn.cluster import KMeans

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
    "E175B": "E175",
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


# %%
def main(idx):
    scale_method = "logcpm"
    count_path = Path.joinpath(
        WORKDIR,
        f"Data/scale_df/{scale_method}-hotspot/{idx}-{scale_method}-hotspot.csv",
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

    selected_genes = []
    with open(
            Path.joinpath(
                WORKDIR,
                f"results/I-gmm/{idx_tp[idx]}-{scale_method}-3.csv",
            )) as f:
        for line in f:
            line = line.strip()
            selected_genes.append(line)
    count_sub_df = count_df.reindex(columns=selected_genes)

    n_spot_clusters = 9

    gene_distmat = sch.distance.pdist(count_sub_df.T, metric="jaccard")
    spot_distmat = sch.distance.pdist(count_sub_df, metric="jaccard")
    Z_gene = sch.linkage(gene_distmat, method="ward")
    gene_result = pd.Series(
        sch.fcluster(Z_gene, t=n_gene_clusters, criterion="maxclust"),
        index=count_sub_df.columns,
    ).sort_values()
    Z_spot = sch.linkage(spot_distmat, method="ward")
    spot_result = pd.Series(
        sch.fcluster(Z_spot, t=n_spot_clusters, criterion="maxclust"),
        index=count_sub_df.index,
    ).sort_values()

    # gene_result = pd.Series(
    # KMeans(n_gene_clusters).fit_predict(count_sub_df.T),
    # index=count_sub_df.columns,
    # ).sort_values() + 1
    # spot_result = pd.Series(
    # KMeans(n_spot_clusters).fit_predict(count_sub_df),
    # index=count_sub_df.index,
    # ).sort_values() + 1

    left, bottom = 0.1, 0.1
    width, height = 0.66, 0.66
    spacing, cluster_width = 0.03, 0.22
    rect_heatmap = [left + spacing, bottom + spacing * 2, width, height]

    fig = plt.figure(figsize=(15, 15))
    ax_heatmap = fig.add_axes(rect_heatmap)

    for i in range(n_gene_clusters):
        rect_cluster = [
            left + width + spacing * (2 + i // 3) + cluster_width * (i // 3),
            bottom + spacing * (3 - i % 3) + cluster_width * (2 - i % 3),
            cluster_width,
            cluster_width,
        ]
        ax_cluster = fig.add_axes(rect_cluster)
        ax_cluster.set_title(f"Cluster {i + 1} hotspot mean", fontsize=10)
        ax_cluster.set_xticks([])
        ax_cluster.set_yticks([])
        ax_cluster.imshow(he_image)
        ax_cluster.axis("off")
        sc = ax_cluster.scatter(
            coor_df["X"],
            coor_df["Y"],
            c=count_sub_df[gene_result[gene_result == i + 1].index].T.mean(),
            cmap="autumn_r",
            vmin=0,
            vmax=0.8,
            s=4,
            alpha=0.7,
        )

    ii = ceil(n_gene_clusters / 3)
    rect_cb = [
        left + width + spacing * (ii + 2) + cluster_width * ii,
        bottom + spacing * 2,
        0.02,
        height,
    ]
    ax_cb = fig.add_axes(rect_cb)
    fig.colorbar(sc, cax=ax_cb)

    ax_heatmap.imshow(
        count_sub_df.reindex(
            columns=gene_result.index,
            index=spot_result.index,
        ),
        cmap="Reds",
        aspect="auto",
    )

    flag = 0
    for i in range(1, n_gene_clusters):
        flag += len(gene_result[gene_result == i])
        ax_heatmap.plot([flag, flag], [0, count_sub_df.shape[0] - 5])
    flag = 0
    for i in range(1, n_spot_clusters):
        flag += len(spot_result[spot_result == i])
        ax_heatmap.plot([0, count_sub_df.shape[1] - 5], [flag, flag])

    ax_heatmap.set_title(f"{idx}")
    ax_heatmap.set_xticks([])
    ax_heatmap.set_yticks([])
    ax_heatmap.set_xlabel("Genes")
    ax_heatmap.set_ylabel("Spots")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/hierarchical/{idx}-{n_gene_clusters}.jpg",
        ),
        bbox_inches="tight",
    )
    plt.close(fig)

    gene_result.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/hierarchical/{idx}-{n_gene_clusters}-genes.csv",
        ))


if __name__ == "__main__":
    for n_gene_clusters in range(6, 15):
        pool = Pool(6)
        pool.map(main, idx_full.keys())
        pool.close()
        pool.join()