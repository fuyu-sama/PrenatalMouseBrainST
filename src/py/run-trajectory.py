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

from arboreto.algo import grnboost2

from pyscenic.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})

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
c = {
    "E135A": "#00FFFF",
    "E135B": "#00FFFF",
    "E155A": "#FFD700",
    "E155B": "#FFD700",
    "E165A": "#0000FF",
    "E165B": "#0000FF",
    "E175A1": "#EE82EE",
    "E175A2": "#EE82EE",
    "E175B": "#EE82EE",
    "P0A1": "#5F9EA0",
    "P0A2": "#5F9EA0",
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
    scale_method = "combat-qq-logcpm"
    cluster_method = "sc3"

# %% read counts
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %%
# TFs from AnimalTFDB3.0
tf_names = []
with open(Path.joinpath(WORKDIR, "Data/SCENIC/Mus_musculus_TF.txt")) as f:
    for line in f:
        tf = line.split("\t")[1]
        if tf == "Symbol":
            continue
        tf_names.append(tf)

db_fnames = [
    Path.joinpath(WORKDIR, "Data/SCENIC", i) for i in [
        "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
        "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
    ]
]
dbs = [RankingDatabase(fname=fname, name=fname) for fname in db_fnames]

#   adjacencies = grnboost2(expression_data=count_full_df, tf_names=tf_names)
#   adjacencies.to_csv(
#       Path.joinpath(WORKDIR, f"results/SCENIC/{scale_method}/adjacencies.csv"))
adjacencies = pd.read_csv(
    Path.joinpath(WORKDIR,
                  f"results/SCENIC/{scale_method}/adjacencies.csv"))

modules = list(modules_from_adjacencies(adjacencies, count_full_df))

def main():
    df = prune2df(
        dbs,
        modules,
        Path.joinpath(WORKDIR, "Data/SCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl"),
    )
    df.to_csv(
        Path.joinpath(WORKDIR, f"results/SCENIC/{scale_method}/motifs.csv"))

if __name__ == "__main__":
    main()
regulons = df2regulons(df)

auc_mtx = aucell(count_full_df, regulons)
auc_mtx.to_csv(
    Path.joinpath(WORKDIR, f"results/SCENIC/{scale_method}/auc.csv"), )
