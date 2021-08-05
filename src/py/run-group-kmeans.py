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
from sklearn.decomposition import PCA, FastICA

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})
Image.MAX_IMAGE_PIXELS = None

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

scale_method = sys.argv[1]

# %% read data
scale_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}-inter.csv",
)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/full-coor.csv")

scale_df = pd.read_csv(scale_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
assert all(scale_df.index == coor_df.index)

he_images = {
    idx: Image.open(Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif"))
    for idx in idx_full
}

# %% pca and kmeans
model_decomposition = PCA(n_components=30)
result_decomposition = model_decomposition.fit_transform(scale_df)
cluster_df = coor_df.copy(deep=True)
for k in range(5, 29):
    model_kmeans = KMeans(n_clusters=k)
    cluster_result = model_kmeans.fit_predict(result_decomposition)
    cluster_df[f"kmeans_{k}_clusters"] = cluster_result
    # draw together
    fig, axes = plt.subplots(3, 4, figsize=(40, 30))
    for idx, ax in zip(idx_full, axes.flatten()):
        draw_df = cluster_df.reindex(
            index=[i for i in cluster_df.index if idx in i])
        ax.imshow(he_images[idx])
        ax.set_title(f"{idx} ncs {k}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.scatter(
            draw_df["X"],
            draw_df["Y"],
            s=16,
            c=draw_df[f"kmeans_{k}_clusters"],
            cmap=ListedColormap(colors[:k]),
            alpha=0.7,
            vmin=0,
            vmax=k,
        )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/group-kmeans/{scale_method}-pca/together/{k}.pdf"
        ))
    plt.close(fig)
    # draw separate
    for idx in idx_full:
        nrows = math.ceil(k / 4)
        cluster_series = cluster_df[f"kmeans_{k}_clusters"].reindex(
            index=[i for i in cluster_df.index if idx in i])
        fig, axes = plt.subplots(nrows, 4, figsize=(40, nrows * 10))
        [i.set_xticks([]) for i in axes.flatten()]
        [i.set_yticks([]) for i in axes.flatten()]
        for n, ax in zip(set(cluster_series), axes.flatten()):
            draw_series = cluster_series[cluster_series == n]
            draw_coor = coor_df.reindex(index=draw_series.index)
            ax.imshow(he_images[idx])
            ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
            ax.set_title(f"{idx} ncs {k} cluster {n}")
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/cluster/group-kmeans/{scale_method}-pca/separate/{idx}-{k}.pdf"
            ))
        plt.close(fig)
cluster_df.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/group-kmeans/{scale_method}-pca/pattern/full-pca-kmeans.csv")
)

# %% ica and kmeans
model_decomposition = FastICA(n_components=30)
result_decomposition = model_decomposition.fit_transform(scale_df)
cluster_df = coor_df.copy(deep=True)
for k in range(5, 29):
    model_kmeans = KMeans(n_clusters=k)
    cluster_result = model_kmeans.fit_predict(result_decomposition)
    cluster_df[f"kmeans_{k}_clusters"] = cluster_result
    # draw together
    fig, axes = plt.subplots(3, 4, figsize=(40, 30))
    for idx, ax in zip(idx_full, axes.flatten()):
        draw_df = cluster_df.reindex(
            index=[i for i in cluster_df.index if idx in i])
        ax.imshow(he_images[idx])
        ax.set_title(f"{idx} ncs {k}")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.scatter(
            draw_df["X"],
            draw_df["Y"],
            s=16,
            c=draw_df[f"kmeans_{k}_clusters"],
            cmap=ListedColormap(colors[:k]),
            alpha=0.7,
            vmin=0,
            vmax=k,
        )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/group-kmeans/{scale_method}-ica/together/{k}.pdf"
        ))
    plt.close(fig)
    # draw separate
    for idx in idx_full:
        nrows = math.ceil(k / 4)
        cluster_series = cluster_df[f"kmeans_{k}_clusters"].reindex(
            index=[i for i in cluster_df.index if idx in i])
        fig, axes = plt.subplots(nrows, 4, figsize=(40, nrows * 10))
        [i.set_xticks([]) for i in axes.flatten()]
        [i.set_yticks([]) for i in axes.flatten()]
        for n, ax in zip(set(cluster_series), axes.flatten()):
            draw_series = cluster_series[cluster_series == n]
            draw_coor = coor_df.reindex(index=draw_series.index)
            ax.imshow(he_images[idx])
            ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
            ax.set_title(f"{idx} ncs {k} cluster {n}")
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/cluster/group-kmeans/{scale_method}-ica/separate/{idx}-{k}.pdf"
            ))
        plt.close(fig)
cluster_df.to_csv(
    Path.joinpath(
        WORKDIR,
        f"results/group-kmeans/{scale_method}-ica/pattern/full-ica-kmeans.csv")
)
