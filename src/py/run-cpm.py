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
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

session_info.show()
plt.rcParams.update({"font.size": 24})
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

# %% read raw count table
idx = sys.argv[1]
he_path = Path.joinpath(
    Path.home(), f"workspace/mouse-brain-full/Data/HE/{idx_full[idx]}.tif"
)
count_path = Path.joinpath(
    Path.home(),
    f"workspace/mouse-brain-full/spaceranger/{idx}/outs/filtered_fature_bc_matrix/{idx}.csv",
)
coor_path = Path.joinpath(
    Path.home(),
    f"workspace/mouse-brain-full/spaceranger/{idx}/outs/spatial/coor-{idx}.csv",
)

he_image = Image.open(he_path)
count_df = pd.read_csv(count_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
count_df.index = [f"{idx}_{i}" for i in count_df.index]
coor_df.index = [f"{idx}_{i}" for i in coor_df.index]
assert all(coor_df.index == count_df.index)

# %% drop low
drop_gene = []
for i in count_df.columns:
    if sum(count_df[i]) <= 0:
        drop_gene.append(i)
count_df.drop(columns=drop_gene, inplace=True)

# %% scale data
scale_df = np.log(count_df.T * 10000 / count_df.T.sum() + 1)
scale_df.to_csv(
    Path.joinpath(
        Path.home(), f"workspace/mouse-brain-full/logcpm/scale_df/{idx}-logcpm.csv"
    )
)
