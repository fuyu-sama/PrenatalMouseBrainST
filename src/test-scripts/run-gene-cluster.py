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
import os
import sys
from pathlib import Path

import libpysal
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
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
    cluster_method = sys.argv[2]
    idx = sys.argv[3]

except IndexError:
    scale_method = "combat-zjq"
    cluster_method = "sc3"
    idx = "E165A"

# %% read data
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}-moran/{idx}-{scale_method}-moran.csv",
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

cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
cluster_df = cluster_df.reindex(
    index=[i for i in cluster_df.index if idx in i])

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

regions_label = {j: i for i, j in enumerate(regions)}
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"{cluster_method}_clusters"]))
    if i not in in_regions
]

draw_cluster = cluster_df[f"{cluster_method}_clusters"].copy()
for c in regions:
    draw_cluster.replace(regions[c], c, inplace=True)
flag = len(regions)
for c in others:
    draw_cluster.replace(c, flag, inplace=True)
for c, i in zip(regions, range(len(regions))):
    c_len = draw_cluster[draw_cluster == c].shape[0]
    draw_cluster.replace(c, regions_label[c], inplace=True)
draw_cluster.sort_values(inplace=True)

# %% pvalue
pvals = SpaGene.spagene.Spagene(
    coordinate_df=coor_df,
    expression_df=count_df,
    genes=count_df.columns,
    permutation=9999,
    cores=20,
)
pvals_cbd = SpaGene.utils.combine_p_values(pvals).sort_values(
    by="adjusted_p_value")
count_sub_df = count_df.reindex(columns=pvals_cbd.index[:1000])

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
ja = 2
ja_path = Path.joinpath(WORKDIR, f"results/5/{idx}/1/ja_{ja}")
ja_set = gene_result[gene_result == ja].index
if not os.path.exists(ja_path):
    os.mkdir(ja_path)
for i in ja_set:
    os.system(f"cp draw_genes/all/{i}.jpg {ja_path}")

# %%
fig, ax = plt.subplots(figsize=(10, 10))
ax.violinplot([
    pvals_cbd.loc[gene_result[gene_result == i].index, "adjusted_p_value"]
    for i in range(1, 1 + n_gene_clusters)
])
fig.savefig(Path.joinpath(WORKDIR, f"results/5/{idx}/jaccard/pvals-adj.jpg"))

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

# %%
points = np.array(coor_df[["X", "Y"]])
weight = libpysal.weights.KNN(points, k=4)
global_moran_df = SpaGene.moran.global_moran(
    gene_lists=count_df.columns,
    expression_df=count_df,
    select_weights=weight,
    transform_moran='r',
    permutation=999,
    cores=20,
)
global_moran_df = global_moran_df.sort_values(by="I_value", ascending=False)

# %%
global_moran_df2 = pd.read_csv(
    Path.joinpath(WORKDIR, "Data/E165A_cpm_knn_gaussian_global_moran_999.csv"),
    index_col=0, header=0,
).sort_values(by="knn4_I_value", ascending=False)

# %%
n_genes = 2000
fig, ax = plt.subplots(figsize=(10, 10))
venn2(
    [set(global_moran_df.index[:n_genes]),
     set(global_moran_df2.index[:n_genes])],
    set_labels=("hotspot", "cpm"),
    ax=ax,
)
ax.set_title(f"Top {n_genes} genes")
fig.savefig(Path.joinpath(WORKDIR, f"results/5/{idx}/{n_genes}.jpg"))
