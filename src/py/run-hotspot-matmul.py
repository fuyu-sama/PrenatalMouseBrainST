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
import sys
from collections import Counter
from pathlib import Path
from PIL import Image

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")

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
    idx = sys.argv[1]
    ncs = int(sys.argv[2])
except IndexError:
    idx = "P0A1"
    ncs = 14

# %% read data
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/logcpm-hotspot-6/{idx}-logcpm-hotspot-6.csv",
)
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

count_df = count_df.astype("int")

assert all(count_df.index == coor_df.index)

# %% read cluster result
spot_cluster_df = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/logcpm-hotspot-6-logcpm-0.99-Ai-0_500",
        f"{idx}-{ncs}/{idx}-spots.csv",
    ),
    index_col=0,
    header=0,
).reindex(index=coor_df.index)

gene_cluster_series = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/logcpm-hotspot-6-logcpm-0.99-Ai-0_500",
        f"{idx}-{ncs}/{idx}-genes.csv",
    ),
    index_col=0,
    header=0,
)["0"]

# %% find neighbor clusters
certain_series = spot_cluster_df[spot_cluster_df["spot_type"] != "uncertain"]
certain_series = certain_series["type_1"]
nbrs = NearestNeighbors(n_neighbors=7).fit(coor_df)
neighbor_clusters_ = {}
for present_cluster in set(spot_cluster_df["type_1"]):
    used_spots = []
    neighbor_clusters_[present_cluster] = []
    present_spots = certain_series[certain_series == present_cluster].index
    distances, indices = nbrs.kneighbors(coor_df.reindex(index=present_spots))
    for record in indices:
        for neighbor_spot in record[1:]:
            if neighbor_spot in used_spots:
                continue
            used_spots.append(neighbor_spot)
            neighbor_cluster = spot_cluster_df.iloc[neighbor_spot, 1]
            if neighbor_cluster == present_cluster:
                continue
            neighbor_clusters_[present_cluster].append(neighbor_cluster)
neighbor_clusters = {}
for i in neighbor_clusters_:
    counter = Counter(neighbor_clusters_[i])
    coverage_ratio = len(certain_series[certain_series == i])
    coverage_ratio = coverage_ratio * 0.05
    neighbor_clusters[i] = [
        j[0] for j in counter.most_common()[:3] if j[1] >= coverage_ratio
    ]
    neighbor_clusters[i].append(i)

# %% draw neighbor clusters
he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = Image.open(he_path)
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(
    coor_df["X"],
    coor_df["Y"],
    s=16,
    c=spot_cluster_df["type_1"],
    cmap='tab20',
)
ax.imshow(he_image)
fig.savefig(Path.joinpath(WORKDIR, f"results/matmul/{idx}-spots.jpg"))

for center_cluster in set(spot_cluster_df["type_1"]):
    draw_clusters = list(neighbor_clusters[center_cluster])
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("off")
    ax.imshow(he_image)
    for i, cluster in enumerate(draw_clusters):
        draw_spots = spot_cluster_df[spot_cluster_df["type_1"] ==
                                     cluster].index
        ax.scatter(
            coor_df.reindex(index=draw_spots)["X"],
            coor_df.reindex(index=draw_spots)["Y"],
            label=cluster,
            s=16,
            c=colors[i],
        )
    draw_spots = spot_cluster_df[spot_cluster_df["type_1"] == center_cluster]
    draw_spots = draw_spots.index
    ax.scatter(
        coor_df.reindex(index=draw_spots)["X"],
        coor_df.reindex(index=draw_spots)["Y"],
        label=center_cluster,
        s=16,
        c="black",
    )
    ax.legend()
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/matmul/{idx}-{center_cluster}-neighbors.jpg",
        ))
    plt.close(fig)

# %%
center_cluster = 10
use_neighbor = True
if use_neighbor:
    selected_spots = [
        j for i in neighbor_clusters[center_cluster]
        for j in spot_cluster_df[spot_cluster_df["type_1"] == i].index
    ]
    selected_genes = [
        j for i in neighbor_clusters[center_cluster]
        for j in gene_cluster_series[gene_cluster_series == i].index
    ]
else:
    selected_spots = spot_cluster_df[spot_cluster_df["type_1"] ==
                                     center_cluster].index
    selected_genes = gene_cluster_series[gene_cluster_series ==
                                         center_cluster].index

count_sub = count_df.reindex(index=selected_spots, columns=selected_genes)
con_matrix = count_sub.T @ count_sub
con_matrix = con_matrix - np.diag(np.diag(con_matrix.to_numpy()))
men_matrix = (1 - count_sub.T) @ count_sub

con_threshold = len(
    spot_cluster_df[spot_cluster_df["type_1"] == center_cluster]) * 0.03
con_matrix_bi = con_matrix.where(con_matrix > con_threshold, 0)
con_matrix_bi = con_matrix_bi.where(con_matrix_bi < 1, 1)
mul_matrix = con_matrix_bi * men_matrix

# %%
gene_series_dict = {}
for gene in mul_matrix.columns:
    gene_series = mul_matrix[gene].sort_values(ascending=False)
    gene_series = gene_series[gene_series > gene_series.mean()]
    gene_series_dict[gene] = gene_series / gene_series.max()

gene_pairs_dict = {}
for gene_1 in gene_series_dict:
    for gene_2 in gene_series_dict[gene_1].index:
        if gene_1 in gene_series_dict[gene_2].index:
            power = (gene_series_dict[gene_1][gene_2] +
                     gene_series_dict[gene_2][gene_1]) / 2
            gene_pair_set = frozenset((gene_1, gene_2))
            if gene_pair_set not in gene_pairs_dict:
                gene_pairs_dict[gene_pair_set] = power

gene_pairs_df = pd.DataFrame(columns=["gene_1", "gene_2", "power"])
for i in gene_pairs_dict:
    temp_series = pd.Series({
        "gene_1": list(i)[0],
        "gene_2": list(i)[1],
        "power": gene_pairs_dict[i],
    })
    gene_pairs_df = pd.concat([gene_pairs_df, temp_series.to_frame().T])
