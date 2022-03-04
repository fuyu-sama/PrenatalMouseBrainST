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

import libpysal
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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
except IndexError:
    scale_method = "cpm"

# %% read counts
count_dict = {}
coor_dict = {}
moran_dict = {}
for idx in idx_full:
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
    count_df = count_df.reindex(index=coor_df.index)

    weights = libpysal.weights.KNN(coor_df, k=8)
    global_moran = SpaGene.moran.global_moran(
        selected_genes=count_df.columns,
        gene_expression_df=count_df,
        weights=weights,
        transform_moran="r",
        permutation=999,
        cores=20,
    )
    global_moran.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/global_moran/{idx}-{scale_method}-8.csv",
        ))
    moran_dict[idx] = global_moran

    count_dict[idx] = count_df
    coor_dict[idx] = coor_df

# %%
for idx in idx_full:
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.hist(
        np.log(moran_dict[idx]["I_value"].dropna().astype(np.float) + 2),
        bins=100,
    )
    ax.set_title(f"{idx} log(I + 2)")
    ax.set_yticks([])
    fig.savefig(Path.joinpath(WORKDIR, f"results/1/{idx}-1.jpg"))

# %%
for idx in idx_full:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    ax1.hist(
        np.log(moran_dict[idx]["I_value"][
            moran_dict[idx]["I_value"] > 0].dropna().astype(np.float)),
        bins=100,
    )
    ax1.set_title(f"{idx} log(I) where I > 0")
    ax1.set_yticks([])
    ax2.hist(
        np.log(-moran_dict[idx]["I_value"][
            moran_dict[idx]["I_value"] < 0].dropna().astype(np.float)),
        bins=100,
    )
    ax2.set_title(f"{idx} log(-I) where I < 0")
    ax2.set_yticks([])
    fig.savefig(Path.joinpath(WORKDIR, f"results/1/{idx}-2.jpg"))