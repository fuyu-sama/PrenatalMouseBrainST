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
import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from PIL import Image
from scipy.cluster import hierarchy as sch

import SpaGene

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
    idx = sys.argv[2]

except IndexError:
    scale_method = "cpm-moran-8"
    idx = "E165A"

# %% read data
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

he_image = Image.open(Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif"))

drop_genes = []
for i in count_df:
    if count_df[i].var() < 1e-6:
        drop_genes.append(i)
count_df.drop(columns=drop_genes, inplace=True)

# %% pvalue
global_moran_df = pd.read_csv(
    Path.joinpath(WORKDIR, "Data/E165A_cpm_knn_gaussian_global_moran_999.csv"),
    index_col=0, header=0,
).sort_values(by="knn8_I_value", ascending=False)
count_sub_df = count_df.reindex(columns=global_moran_df.index[:1000])

# %% calculate distmat
gene_distmat = sch.distance.pdist(count_sub_df.T, metric="jaccard")
spot_distmat = sch.distance.pdist(count_sub_df, metric="jaccard")

# %% cluster
n_gene_clusters = 12
n_spot_clusters = 12

Z = sch.linkage(gene_distmat, method="ward")
gene_result = pd.Series(
    sch.fcluster(Z, t=n_gene_clusters, criterion="maxclust"),
    index=count_sub_df.columns,
).sort_values()

Z = sch.linkage(spot_distmat, method="ward")
spot_result = pd.Series(
    sch.fcluster(Z, t=n_spot_clusters, criterion="maxclust"),
    index=count_sub_df.index,
).sort_values()

# %% draw
fig, ax = plt.subplots(figsize=(10, 10))
ax.imshow(
    count_sub_df.reindex(columns=gene_result.index, index=spot_result.index),
    cmap="Reds",
    aspect="auto",
)
ax.set_title("jaccard")
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("Genes")
ax.set_ylabel("Spots")
flag = 0
for i in range(1, n_gene_clusters):
    flag += len(gene_result[gene_result == i])
    ax.plot([flag, flag], [0, count_df.shape[0] - 5])
fig.savefig(Path.joinpath(WORKDIR, f"results/5/{idx}/jaccard/jaccard-0.jpg"),
            bbox_inches="tight")

for i in range(1, n_gene_clusters + 1):
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.imshow(he_image)
    sc = ax.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=count_df[gene_result[gene_result == i].index].T.mean(),
        cmap="autumn_r",
        s=16,
        alpha=0.7,
    )
    ax.set_title(f"Cluster {i} jaccard genes")
    fig.colorbar(sc, ax=ax)
    fig.savefig(
        Path.joinpath(WORKDIR, f"results/5/{idx}/jaccard/jaccard-{i}.jpg"),
        bbox_inches="tight",
    )

# %%
for gene in ["Gbx2", "Zbtb18", "Calb2", "Satb2", "Tbr1", "Sox5", "Sox2"]:
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=count_df[gene],
        cmap="Reds",
        vmin=0,
        vmax=1,
    )
    ax.imshow(he_image)
    ax.set_title(gene)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(Path.joinpath(WORKDIR, f"results/5/{idx}/{gene}.jpg"))
