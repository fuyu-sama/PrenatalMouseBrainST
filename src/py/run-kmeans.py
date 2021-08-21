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
import math
import sys
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image
from sklearn.cluster import KMeans
from sklearn.decomposition import FastICA, PCA

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})
Image.MAX_IMAGE_PIXELS = None
random_state = 42

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
    "P0B": "V10M17-100-P0B",
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

idx = sys.argv[1]
scale_method = sys.argv[2]
raw_scale = int(sys.argv[3])
pca = int(sys.argv[4])
ica = int(sys.argv[5])

# %% read data
he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
scale_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")

he_image = Image.open(he_path)
scale_df = pd.read_csv(scale_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
assert all(scale_df.index == coor_df.index)

# %% kmeans
if raw_scale:
    cluster_df = pd.DataFrame(index=scale_df.index)
    fig, ax = plt.subplots(6, 4, figsize=(40, 60))
    for k, j in zip(range(5, 29), ax.flatten()):
        model_kmeans = KMeans(n_clusters=k, random_state=random_state)
        cluster_result = model_kmeans.fit_predict(scale_df)
        cluster_df[f"kmeans_{k}_clusters"] = cluster_result
        j.axis("off")
        j.imshow(he_image)
        j.set_title(f"{idx} ncs {k}")
        j.scatter(
            coor_df["X"],
            coor_df["Y"],
            s=16,
            c=cluster_result,
            cmap=ListedColormap(colors[:k]),
            alpha=0.7,
        )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-kmeans/together/{idx}-kmeans.jpg"
        ))
    plt.close(fig)
    cluster_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-kmeans/pattern/{idx}-kmeans.csv"
        ))

# %% draw separate
if raw_scale:
    for k in range(5, 29):
        cluster_series = cluster_df[f"kmeans_{k}_clusters"]
        nrows = math.ceil(k / 5)
        fig, axes = plt.subplots(nrows, 5, figsize=(50, nrows * 10))
        [ax.axis("off") for ax in axes.flatten()]
        for n, ax in zip(set(cluster_series), axes.flatten()):
            draw_series = cluster_series[cluster_series == n]
            draw_coor = coor_df.reindex(index=draw_series.index)
            ax.imshow(he_image)
            ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
            ax.set_title(f"{idx} ncs {k} cluster {n}")
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster/{scale_method}-kmeans/separate/{idx}-{k}-kmeans.jpg"
                ))
        plt.close(fig)

# %% pca and kmeans
if pca:
    model_decomposition = PCA(n_components=30, random_state=random_state)
    result_decomposition = model_decomposition.fit_transform(scale_df)
    cluster_df = pd.DataFrame(index=scale_df.index)
    fig, ax = plt.subplots(6, 4, figsize=(40, 60))
    for k, j in zip(range(5, 29), ax.flatten()):
        model_kmeans = KMeans(n_clusters=k, random_state=random_state)
        cluster_result = model_kmeans.fit_predict(result_decomposition)
        cluster_df[f"pca-kmeans_{k}_clusters"] = cluster_result
        j.axis("off")
        j.imshow(he_image)
        j.set_title(f"{idx} ncs {k}")
        j.scatter(
            coor_df["X"],
            coor_df["Y"],
            s=16,
            c=cluster_result,
            cmap=ListedColormap(colors[:k]),
            alpha=0.7,
        )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-pca-kmeans/together/{idx}-pca-kmeans.jpg"
        ))
    plt.close(fig)
    cluster_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-pca-kmeans/pattern/{idx}-pca-kmeans.csv"
        ))

# %% draw separate
if pca:
    for k in range(5, 29):
        cluster_series = cluster_df[f"pca-kmeans_{k}_clusters"]
        nrows = math.ceil(k / 5)
        fig, axes = plt.subplots(nrows, 5, figsize=(50, nrows * 10))
        [ax.axis("off") for ax in axes.flatten()]
        for n, ax in zip(set(cluster_series), axes.flatten()):
            draw_series = cluster_series[cluster_series == n]
            draw_coor = coor_df.reindex(index=draw_series.index)
            ax.imshow(he_image)
            ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
            ax.set_title(f"{idx} ncs {k} cluster {n}")
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster/{scale_method}-pca-kmeans/separate/{idx}-{k}-pca-kmeans.jpg"
                ))
        plt.close(fig)


# %% ica and kmeans
if ica:
    model_decomposition = FastICA(n_components=30, random_state=random_state)
    result_decomposition = model_decomposition.fit_transform(scale_df)
    cluster_df = pd.DataFrame(index=scale_df.index)
    fig, ax = plt.subplots(6, 4, figsize=(40, 60))
    for k, j in zip(range(5, 29), ax.flatten()):
        model_kmeans = KMeans(n_clusters=k, random_state=random_state)
        cluster_result = model_kmeans.fit_predict(result_decomposition)
        cluster_df[f"ica-kmeans_{k}_clusters"] = cluster_result
        j.axis("off")
        j.imshow(he_image)
        j.set_title(f"{idx} ncs {k}")
        j.scatter(
            coor_df["X"],
            coor_df["Y"],
            s=16,
            c=cluster_result,
            cmap=ListedColormap(colors[:k]),
            alpha=0.7,
        )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-ica-kmeans/together/{idx}-ica-kmeans.jpg"
        ))
    plt.close(fig)
    cluster_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-ica-kmeans/pattern/{idx}-ica-kmeans.csv"
        ))

# %% draw separate
if ica:
    for k in range(5, 29):
        cluster_series = cluster_df[f"ica-kmeans_{k}_clusters"]
        nrows = math.ceil(k / 5)
        fig, axes = plt.subplots(nrows, 5, figsize=(50, nrows * 10))
        [ax.axis("off") for ax in axes.flatten()]
        for n, ax in zip(set(cluster_series), axes.flatten()):
            draw_series = cluster_series[cluster_series == n]
            draw_coor = coor_df.reindex(index=draw_series.index)
            ax.imshow(he_image)
            ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
            ax.set_title(f"{idx} ncs {k} cluster {n}")
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster/{scale_method}-ica-kmeans/separate/{idx}-{k}-ica-kmeans.jpg"
                ))
        plt.close(fig)
