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
from sklearn.decomposition import PCA

plt.rcParams.update({"font.size": 24})
Image.MAX_IMAGE_PIXELS = None
random_state = 42
session_info.show()

idx = sys.argv[1]

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
    "#FAEBD7",
    "#00FFFF",
    "#FFD700",
    "#0000FF",
    "#FF8C00",
    "#EE82EE",
    "#9ACD32",
    "#5F9EA0",
    "#7FFF00",
    "#7FFFD4",
    "#6495ED",
    "#008B8B",
    "#B8860B",
    "#C0C0C0",
    "#000080",
    "#D8BFD8",
    "#00CED1",
    "#9400D3",
    "#8E804B",
    "#0089A7",
    "#CB1B45",
    "#FFB6C1",
    "#00FF00",
    "#800000",
    "#376B6D",
    "#D8BFD8",
    "#F5F5F5",
    "#D2691E",
]

# %% read data
he_path = Path.joinpath(
    Path.home(), f"workspace/mouse-brain-full/Data/HE/{idx_full[idx]}.tif"
)
count_path = Path.joinpath(
    Path.home(), f"workspace/mouse-brain-full/logcpm/scale_df/{idx}-logcpm-inter.csv"
)
coor_path = Path.joinpath(
    Path.home(),
    f"workspace/mouse-brain-full/spaceranger/{idx}/outs/spatial/coor-{idx}.csv",
)

he_image = Image.open(he_path)
count_df = pd.read_csv(count_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
coor_df.index = [f"{idx}_{i}" for i in coor_df.index]
# count_df = count_df.reindex(index=coor_df.index)
assert all(count_df.index == coor_df.index)

# %% clustering and draw
cluster_df = pd.DataFrame(index=count_df.index)
ncs_group = range(5, 29)
fig, axes = plt.subplots(
    6,
    int(len(ncs_group) / 6),
    figsize=(int(len(ncs_group) / 6) * 10, 60),
)
for ncs, ax in zip(ncs_group, axes.flatten()):
    [axis.set_visible(False) for axis in ax.spines.values()]
    model_pca = PCA(n_components=30, random_state=random_state)
    model_kmeans = KMeans(n_clusters=ncs, random_state=random_state)
    result_pca = model_pca.fit_transform(count_df)
    cluster_array = model_kmeans.fit_predict(result_pca)
    cluster_df[f"kmeans_{ncs}_clusters"] = cluster_array
    ax.imshow(he_image)
    ax.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=cluster_df[f"kmeans_{ncs}_clusters"],
        cmap=ListedColormap(colors[:ncs]),
        alpha=0.7,
        s=16,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{idx} {ncs}")
fig.savefig(
    Path.joinpath(
        Path.home(), f"workspace/mouse-brain-full/logcpm/kmeans/cluster/{idx}.jpg"
    )
)
plt.close(fig)
cluster_df.to_csv(
    Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/logcpm/kmeans/pattern/{idx}_KMeans.csv",
    )
)

# %% draw separate
for ncs in ncs_group:
    cluster_series = cluster_df[f"kmeans_{ncs}_clusters"]
    nrows = math.ceil(ncs / 5)
    fig, axes = plt.subplots(nrows, 5, figsize=(50, nrows * 10))
    [ax.set_xticks([]) for ax in axes.flatten()]
    [ax.set_yticks([]) for ax in axes.flatten()]
    for n, ax in zip(set(cluster_series), axes.flatten()):
        draw_series = cluster_series[cluster_series == n]
        draw_coor = coor_df.reindex(index=draw_series.index)
        [axis.set_visible(False) for axis in ax.spines.values()]
        ax.imshow(he_image)
        ax.scatter(draw_coor["X"], draw_coor["Y"], s=16)
        ax.set_title(f"{idx} KMeans cluster {n}")
    fig.savefig(
        Path.joinpath(
            Path.home(),
            f"workspace/mouse-brain-full/logcpm/kmeans/cluster/{idx}-{ncs}.pdf",
        )
    )
    plt.close(fig)
