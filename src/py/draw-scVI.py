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

# %% envoronment config
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import session_info
from scipy import stats
from matplotlib import pyplot as plt

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()

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
raw_df = pd.read_csv(
    Path.joinpath(WORKDIR, f"Data/scale_df/raw_count/{idx}-raw.csv"),
    index_col=0,
    header=0,
)
scvi_df = pd.read_csv(
    Path.joinpath(WORKDIR, f"Data/scale_df/scvi/{idx}-scvi.csv"),
    index_col=0,
    header=0,
)

pearson_df = pd.DataFrame(index=scvi_df.index, columns=scvi_df.index)
for i in range(pearson_df.shape[0]):
    for j in range(i, pearson_df.shape[0]):
        pearson_df.iloc[i, j] = np.corrcoef(raw_df.iloc[i], scvi_df.iloc[j])[1, 0]

spearman_df = pd.DataFrame(index=scvi_df.index, columns=scvi_df.index)
for i in range(spearman_df.shape[0]):
    for j in range(i, spearman_df.shape[0]):
        spearman_df.iloc[i, j] = stats.spearmanr(raw_df.iloc[i], scvi_df.iloc[j])[1]

fig, ax = plt.subplots(1, 2, figsize=(10, 20))
ax[0].imshow(pearson_df, cmap="bwr")
ax[0].set_title(f"{idx} Pearson")
ax[1].imshow(spearman_df, cmap="bwr")
ax[1].set_title(f"{idx} Spearman")
[i.set_xticks([]) for i in ax]
[i.set_yticks([])for y in ax]
fig.savefig(Path.joinpath(WORKDIR, f"Data/scale_df/scvi/{idx}.jpg"))
