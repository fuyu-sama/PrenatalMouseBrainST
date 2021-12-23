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
from pathlib import Path

import libpysal
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from umap import UMAP

from SpaGene.weighted_correlation import Moran_Process

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
except IndexError:
    scale_method = "combat-zjq"
    cluster_method = "sc3"

# %% read counts
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
    points = np.array(coor_df[["X", "Y"]])
    weight = libpysal.weights.KNN(points, k=4)
    moran_df = Moran_Process(
        selected_genes=count_df.columns,
        gene_expression_df=count_df,
        weights_moran=weight,
        transform_moran='r',
        permutation=999,
        cores=20,
    )
    moran_dict[idx] = moran_df

    dim_results = UMAP().fit_transform(moran_df)
    dim_results = pd.DataFrame(dim_results, index=moran_df.index)

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title(idx)
    ax.scatter(
        dim_results[0],
        dim_results[1],
        s=4,
    )
    fig.savefig(Path.joinpath(WORKDIR, f"results/3/{idx}.jpg"))
    plt.close(fig)