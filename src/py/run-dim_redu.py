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
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from PIL import Image
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

try:
    scale_method = sys.argv[1]
except IndexError:
    scale_method = "combat-Ai-0_500"

# %% read cluster result
region_path = Path.joinpath(
    WORKDIR, f"results/gene-cluster/{scale_method}/regions.json")
with open(region_path) as f:
    regions = json.load(f)
with open(Path.joinpath(WORKDIR, f"results/cluster/colors.json")) as f:
    colors = json.load(f)

for idx in regions:
    ncs = len(regions[idx])
    cluster_path = Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/{scale_method}/{idx}-{ncs}/{idx}-spots.csv",
    )
    cluster_series = pd.read_csv(
        cluster_path,
        index_col=0,
        header=0,
    )[f"1_type"].astype(str)
    for i in regions[idx]:
        if regions[idx][i] in colors:
            flag = colors[regions[idx][i]]
        else:
            flag = regions[idx][i]
        cluster_series.replace(i, flag, inplace=True)

    count_path = Path.joinpath(
        WORKDIR,
        f"Data/scale_df/combat/{idx}-combat.csv",
    )
    count_df = pd.read_csv(
        count_path,
        index_col=0,
        header=0,
    ).T

    coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)

    assert all(count_df.index == coor_df.index)

    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    he_image = Image.open(he_path)

    write_dir = Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}",
    )
    if not write_dir.exists():
        write_dir.mkdir()

    # tSNE
    result_tsne = TSNE().fit_transform(
        PCA(n_components=40).fit_transform(count_df))
    result_df = pd.DataFrame(
        result_tsne,
        index=count_df.index,
        columns=["X", "Y"],
    )
    result_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-tsne.csv",
        ))

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.scatter(
        result_df["X"],
        result_df["Y"],
        c=cluster_series,
        s=8,
    )
    ax.set_title(f"tSNE {scale_method}")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-tsne-1.jpg",
        ))
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(23, 10))
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.scatter(
        result_df["X"],
        result_df["Y"],
        c=cluster_series,
        s=8,
    )
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=cluster_series,
        s=16,
        alpha=0.7,
    )
    ax2.imshow(he_image)
    fig.suptitle(f"tSNE {scale_method}")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-tsne-2.jpg",
        ))
    plt.close(fig)

    # UMAP
    result_umap = UMAP().fit_transform(
        PCA(n_components=40).fit_transform(count_df))
    result_df = pd.DataFrame(
        result_umap,
        index=count_df.index,
        columns=["X", "Y"],
    )
    result_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-umap.csv",
        ))

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.scatter(
        result_df["X"],
        result_df["Y"],
        c=cluster_series,
        s=8,
    )
    ax.set_title(f"UMAP {scale_method}")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-umap-1.jpg",
        ))
    plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(23, 10))
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.scatter(
        result_df["X"],
        result_df["Y"],
        c=cluster_series,
        s=8,
    )
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=cluster_series,
        s=16,
        alpha=0.7,
    )
    ax2.imshow(he_image)
    fig.suptitle(f"UMAP {scale_method}")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/dimension_reduction/{scale_method}/{idx}-umap-2.jpg",
        ))
    plt.close(fig)

    he_image.close()
