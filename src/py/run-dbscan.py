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
import random
import sys
from pathlib import Path

import libpysal as ps
import pandas as pd
import session_info
from matplotlib import pyplot as plt
from scipy.spatial.distance import pdist
from pointpats import PointPattern
import pointpats.quadrat_statistics as qs


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

try:
    idx = sys.argv[1]
except IndexError:
    idx = "E165A"

# %% read data
count_path = Path.joinpath(WORKDIR, f"Data/scale_df/raw/{idx}-raw.csv")
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

# %% genarate pseudo points
f = open("results/dbscan.csv", "w")
f.write("gene,chi2,df,chi2_pvalue\n")
for gene in count_df.columns:
    meta_dist = min(pdist(coor_df))
    sigma = meta_dist / 3
    pseudo_point = pd.DataFrame(columns=["X", "Y"])
    for i in coor_df.index:
        for j in range(count_df.loc[i, gene] * 10 + 1):
            pseudo_point.loc[f"{i}-{j}", "X"] = random.gauss(
                coor_df.loc[i, "X"],
                sigma,
            )
            pseudo_point.loc[f"{i}-{j}", "Y"] = random.gauss(
                coor_df.loc[i, "Y"],
                sigma,
            )

    # fig, ax = plt.subplots(1, 2, figsize=(60, 30))
    # [i.set_xlim(min(coor_df["X"]) - 200, max(coor_df["X"] + 200)) for i in ax]
    # [i.set_ylim(min(coor_df["Y"]) - 200, max(coor_df["Y"] + 200)) for i in ax]
    # ax[0].scatter(pseudo_point["X"], pseudo_point["Y"], s=1)
    # ax[1].scatter(coor_df["X"], coor_df["Y"], c=count_df[gene], cmap="Reds", vmin=-1)
    # fig.savefig("test.jpg")

    pp = PointPattern(pseudo_point)
    q_h = qs.QStatistic(pp, shape="hexagon")
    f.write(f"{gene},{q_h.chi2},{q_h.df},{q_h.chi2_pvalue}\n")

f.close()
