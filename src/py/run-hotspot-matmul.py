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
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture
from scipy.stats import norm

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
plt.rcParams.update({"font.size": 18})

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
    coverage_ratio = coverage_ratio * 0.03
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
fig.savefig(Path.joinpath(WORKDIR, f"results/matmul.2/{idx}-spots.svg"))

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
            f"results/matmul.2/{idx}-{center_cluster}-neighbors.svg",
        ))
    plt.close(fig)

# %%
gene_pairs = {}
for center_cluster in range(1, 1 + ncs):
    selected_spots = [
        j for i in neighbor_clusters[center_cluster]
        for j in spot_cluster_df[spot_cluster_df["type_1"] == i].index
    ]
    selected_genes = [
        j for i in neighbor_clusters[center_cluster]
        for j in gene_cluster_series[gene_cluster_series == i].index
    ]

    count_sub = count_df.reindex(index=selected_spots, columns=selected_genes)
    con_matrix = count_sub.T @ count_sub
    n_hotspots = np.diag(con_matrix.to_numpy())
    con_matrix = con_matrix - np.diag(n_hotspots)
    men_matrix = (1 - count_sub.T) @ count_sub

    con_ratio = con_matrix / len(selected_spots)
    men_ratio = men_matrix / len(selected_spots)

    gene_pairs_df = pd.DataFrame(columns=[
        "SVG_cluster", "gene_1", "gene_2", "colocalization_score",
        "exclusive_score"
    ])
    line = 0
    for gene_1 in con_ratio.columns:
        for gene_2 in con_ratio.index:
            if gene_1 == gene_2:
                continue
            write_series = pd.Series(
                {
                    "SVG_cluster": center_cluster,
                    "gene_1": gene_1,
                    "gene_2": gene_2,
                    "colocalization_score": con_ratio[gene_1][gene_2],
                    "exclusive_score": men_ratio[gene_1][gene_2],
                },
                name=line,
            ).to_frame().T
            line += 1
            gene_pairs_df = pd.concat([gene_pairs_df, write_series])
    gene_pairs[center_cluster] = gene_pairs_df

# %%
gmm_result = {}
gmm_model = {}
n_components = 3
for center_cluster in range(1, 1 + ncs):
    gmm_result[center_cluster] = {}
    co_df = gene_pairs[center_cluster][
        gene_pairs[center_cluster]["colocalization_score"] > 0].drop(
            columns=["exclusive_score"])
    ex_df = gene_pairs[center_cluster][
        gene_pairs[center_cluster]["exclusive_score"] > 0].drop(
            columns=["colocalization_score"])
    co_em_df = co_df.copy()
    ex_em_df = ex_df.copy()
    co_bayesian_df = co_df.copy()
    ex_bayesian_df = ex_df.copy()
    fig, (
        (ax_co_1, ax_ex_1),
        (ax_co_2, ax_ex_2),
    ) = plt.subplots(2, 2, figsize=(20, 20))

    # colocalization Gaussian mixture with EM
    ax_co_1.hist(
        co_df["colocalization_score"],
        bins=100,
        density=True,
        color="#f0a1a8",
    )
    ax_co_1.set_title(
        f"Cluster {center_cluster} colocalization score distribution, EM")
    gmm_co = GaussianMixture(
        n_components=n_components,
        max_iter=1000,
        n_init=10,
    )
    co_em_df["co_em"] = gmm_co.fit_predict(
        co_df["colocalization_score"].to_numpy().reshape(-1, 1))
    # co_em_df = co_em_df[co_em_df["co_em"] == gmm_co.means_.argmax()]

    x = np.linspace(
        co_df["colocalization_score"].min(),
        co_df["colocalization_score"].max(),
        5000,
    )
    ax_co_1.plot(
        x,
        np.exp(gmm_co.score_samples(x.reshape(-1, 1))),
        lw=4,
        label="GMM, EM",
    )
    for i in range(n_components):
        ax_co_1.plot(
            x,
            norm.pdf(
                x,
                gmm_co.means_[i, 0],
                gmm_co.covariances_[i, 0]**(1 / 2),
            ) * gmm_co.weights_[i],
            lw=4,
            ls="--",
            label=f"Gaussian {i}, weight {gmm_co.weights_[i]:.3f}",
        )
    ax_co_1.legend()

    # exclusive with EM
    ax_ex_1.hist(
        ex_df["exclusive_score"],
        bins=100,
        density=True,
        color="#f0a1a8",
    )
    ax_ex_1.set_title(
        f"Cluster {center_cluster} exclusive score distribution, EM")
    gmm_ex = GaussianMixture(
        n_components=n_components,
        max_iter=1000,
        n_init=10,
    )
    ex_em_df["ex_em"] = gmm_ex.fit_predict(
        ex_df["exclusive_score"].to_numpy().reshape(-1, 1))
    # ex_em_df = ex_em_df[ex_em_df["ex_em"] == gmm_ex.means_.argmax()]

    x = np.linspace(
        ex_df["exclusive_score"].min(),
        ex_df["exclusive_score"].max(),
        5000,
    )
    ax_ex_1.plot(
        x,
        np.exp(gmm_ex.score_samples(x.reshape(-1, 1))),
        lw=4,
        label="GMM, EM",
    )
    for i in range(n_components):
        ax_ex_1.plot(
            x,
            norm.pdf(
                x,
                gmm_ex.means_[i, 0],
                gmm_ex.covariances_[i, 0]**(1 / 2),
            ) * gmm_ex.weights_[i],
            lw=4,
            ls="--",
            label=f"Gaussian {i}, weight {gmm_ex.weights_[i]:.3f}",
        )
    ax_ex_1.legend()

    # colocalizaton Bayesian
    ax_co_2.hist(
        co_df["colocalization_score"],
        bins=100,
        density=True,
        color="#f0a1a8",
    )
    ax_co_2.set_title(
        f"Cluster {center_cluster} colocalization score distribution, Bayesian"
    )
    gmm_co_b = BayesianGaussianMixture(
        n_components=n_components,
        max_iter=1000,
        n_init=10,
    )
    co_bayesian_df["co_bayesian"] = gmm_co_b.fit_predict(
        co_df["colocalization_score"].to_numpy().reshape(-1, 1))
    # co_bayesian_df = co_bayesian_df[co_bayesian_df["co_bayesian"] ==
    # gmm_co_b.means_.argmax()]

    x = np.linspace(
        co_df["colocalization_score"].min(),
        co_df["colocalization_score"].max(),
        5000,
    )
    ax_co_2.plot(
        x,
        np.exp(gmm_co_b.score_samples(x.reshape(-1, 1))),
        lw=4,
        label="GMM, Bayesian",
    )
    for i in range(n_components):
        ax_co_2.plot(
            x,
            norm.pdf(
                x,
                gmm_co_b.means_[i, 0],
                gmm_co_b.covariances_[i, 0]**(1 / 2),
            ) * gmm_co_b.weights_[i],
            lw=4,
            ls="--",
            label=f"Gaussian {i}, weight {gmm_co_b.weights_[i]:.3f}",
        )
    ax_co_2.legend()

    # exclusive Bayesian
    ax_ex_2.hist(
        ex_df["exclusive_score"],
        bins=100,
        density=True,
        color="#f0a1a8",
    )
    ax_ex_2.set_title(
        f"Cluster {center_cluster} exclusive score distribution, Bayesian")
    gmm_ex_b = BayesianGaussianMixture(
        n_components=n_components,
        max_iter=1000,
        n_init=10,
    )
    ex_bayesian_df["ex_bayesian"] = gmm_ex_b.fit_predict(
        ex_df["exclusive_score"].to_numpy().reshape(-1, 1))
    # ex_bayesian_df = ex_bayesian_df[ex_bayesian_df["ex_bayesian"] ==
    # gmm_ex_b.means_.argmax()]

    x = np.linspace(
        ex_df["exclusive_score"].min(),
        ex_df["exclusive_score"].max(),
        5000,
    )
    ax_ex_2.plot(
        x,
        np.exp(gmm_ex_b.score_samples(x.reshape(-1, 1))),
        lw=4,
        label="GMM, Bayesian",
    )
    for i in range(n_components):
        ax_ex_2.plot(
            x,
            norm.pdf(
                x,
                gmm_ex_b.means_[i, 0],
                gmm_ex_b.covariances_[i, 0]**(1 / 2),
            ) * gmm_ex_b.weights_[i],
            lw=4,
            ls="--",
            label=f"Gaussian {i}, weight {gmm_ex_b.weights_[i]:.3f}",
        )
    ax_ex_2.legend()

    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/matmul.2/{idx}-{center_cluster}-dist.svg",
        ))
    plt.close(fig)

    gmm_result[center_cluster] = [
        [co_em_df, ex_em_df],
        [co_bayesian_df, ex_bayesian_df],
    ]
    gmm_model[center_cluster] = [
        [gmm_co, gmm_ex],
        [gmm_co_b, gmm_ex_b],
    ]

# %%
for center_cluster in range(1, 1 + ncs):
    gmm_result[center_cluster][1][0].to_csv(
        Path.joinpath(WORKDIR,
                      f"results/matmul.2/{idx}-{center_cluster}-co.csv"))
    gmm_result[center_cluster][1][1].to_csv(
        Path.joinpath(WORKDIR,
                      f"results/matmul.2/{idx}-{center_cluster}-ex.csv"))
