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

# %% read data
idx = sys.argv[1]
scale_method = sys.argv[2]
he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
scale_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/{idx}-{scale_method}-inter.csv",
)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")

he_image = Image.open(he_path)
scale_df = pd.read_csv(scale_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
assert all(scale_df.index == coor_df.index)

# %% pca and kmeans
model_decomposition = PCA(n_components=30)
result_decomposition = model_decomposition.fit_transform(scale_df)
fig, ax = plt.subplots(6, 4, figsize=(40, 60))
for k, j in zip(range(5, 29), ax.flatten()):
    model_kmeans = KMeans(n_clusters=k)
    cluster_result = model_kmeans.fit_predict(result_decomposition)
    j.imshow(he_image)
    j.set_title(f"{idx} ncs {k}")
    j.set_xticks([])
    j.set_yticks([])
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
        f"results/cluster/{scale_method}-pca-kmeans/together/{idx}-PCA-KMeans.pdf"
    ))
plt.close(fig)

# %% ica and kmeans
model_decomposition = FastICA(n_components=30)
result_decomposition = model_decomposition.fit_transform(scale_df)
fig, ax = plt.subplots(6, 4, figsize=(40, 60))
for k, j in zip(range(5, 29), ax.flatten()):
    model_kmeans = KMeans(n_clusters=k)
    cluster_result = model_kmeans.fit_predict(result_decomposition)
    j.imshow(he_image)
    j.set_title(f"{idx} ncs {k}")
    j.set_xticks([])
    j.set_yticks([])
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
        f"results/cluster/{scale_method}-ica-kmeans/together/{idx}-ICA-KMeans.pdf"
    ))
plt.close(fig)
