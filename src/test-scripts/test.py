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
import os
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from PIL import Image
from scipy.cluster import hierarchy as sch

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

# %% read data
for idx in idx_full:
    count_full_df = pd.DataFrame()
    count_path = Path.joinpath(
        WORKDIR,
        f"spaceranger/{idx}/outs/filtered_feature_bc_matrix/{idx}.csv")
    coor_path = Path.joinpath(
        WORKDIR,
        f"spaceranger/{idx}/outs/spatial/coor-{idx}.csv",
    )

    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    count_df.index = [f"{idx}_{i}" for i in count_df.index]
    coor_df.index = [f"{idx}_{i}" for i in coor_df.index]
    assert all(coor_df.index == count_df.index)

    temp_df = count_df.where(count_df < 1, 1)
    temp_series = temp_df.sum() / temp_df.shape[0]
    drop_genes = []
    for i in temp_series.index:
        if temp_series[i] > 0.99 or temp_series[i] < 0.01:
            drop_genes.append(i)
    count_df.drop(columns=drop_genes, inplace=True)

    drop_genes = []
    for i in count_df.columns:
        if sum(count_df[i]) <= 10:
            drop_genes.append(i)
    count_df.drop(columns=drop_genes, inplace=True)

    count_full_df = pd.concat([count_full_df, count_df], axis=0, join="outer")

count_full_df = count_full_df.dropna(axis="columns")