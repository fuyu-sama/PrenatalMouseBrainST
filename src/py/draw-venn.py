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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
from sklearn.linear_model import LinearRegression

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
tp_full = {
    "E135": ["E135A", "E135B"],
    "E155": ["E155A", "E155B"],
    "E165": ["E165A", "E165B"],
    "E175": ["E175A1", "E175A2"],
    # "E175": ["E175A1", "E175A2", "E175B"],
    "P0": ["P0A1", "P0A2"],
}
idx_tp = {
    "E135A": "E135",
    "E135B": "E135",
    "E155A": "E155",
    "E155B": "E155",
    "E165A": "E165",
    "E165B": "E165",
    "E175A1": "E175",
    "E175A2": "E175",
    # "E175B": "E175",
    # "P0B": "P0",
    "P0A1": "P0",
    "P0A2": "P0",
}
tp_name = {
    "E135": "E13.5",
    "E155": "E15.5",
    "E165": "E16.5",
    "E175": "E17.5",
    "P0": "P0"
}
colors = [
    "#FAEBD7", "#00FFFF", "#FFD700", "#0000FF", "#FF8C00", "#EE82EE",
    "#9ACD32", "#5F9EA0", "#7FFF00", "#7FFFD4", "#6495ED", "#008B8B",
    "#B8860B", "#C0C0C0", "#000080", "#D8BFD8", "#00CED1", "#9400D3",
    "#8E804B", "#0089A7", "#CB1B45", "#FFB6C1", "#00FF00", "#800000",
    "#376B6D", "#D8BFD8", "#F5F5F5", "#D2691E"
]

try:
    tp = sys.argv[1]
except IndexError:
    tp = "P0"

# %%
moran_list = []
for idx in tp_full[tp]:
    moran_list.append(
        pd.read_csv(
            Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
            index_col=0,
            header=0,
        ).sort_values(by="Ai", ascending=False))

# %%
len_inter = []
venn_function = venn2 if len(moran_list) == 2 else venn3
for i in range(0, 5000, 500):
    fig, ax = plt.subplots(figsize=(10, 10))
    set_list = [set(j.iloc[0:i + 500, 0].index) for j in moran_list]
    venn_function(set_list, tp_full[tp], ax=ax)
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/Ai/venn/{tp}-1_{i + 500}.jpg",
        ))
    plt.close(fig)

    for j, k in enumerate(set_list):
        if j == 0:
            inters = k
        else:
            inters = inters & k
    len_inter.append(len(inters) / (i + 500))

fig, ax = plt.subplots(figsize=(13, 10))
ax.plot(range(len(len_inter)), len_inter, "D-")
ax.set_xticks(range(10))
ax.set_xticklabels(
    [f"1 - {i + 500}" for i in range(0, 5000, 500)],
    rotation=45,
)
ax.set_ylim([0, 0.9])
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/Ai/venn/{tp}-plot.jpg",
    ),
    bbox_inches="tight",
)

# %%
fig, ax = plt.subplots(figsize=(10, 10))
x = moran_list[0]["Ai"].to_numpy()
y = moran_list[1]["Ai"].reindex(index=moran_list[0].index).fillna(0).to_numpy()
ax.scatter(x, y, s=4)
reg = LinearRegression().fit(x.reshape(-1, 1), y.reshape(-1, 1))
ax.plot(
    [i / 100 for i in range(100)],
    [(i * reg.coef_[0, 0] + reg.intercept_[0]) / 100 for i in range(100)],
    "r--",
    linewidth=2,
)
pearson = np.corrcoef(
    moran_list[0]["Ai"],
    moran_list[1]["Ai"].reindex(index=moran_list[0].index).fillna(0),
)[0, 1]
reg_func = f"y = {reg.coef_[0, 0]: .4f} * x + {reg.intercept_[0]: .4f}"
ax.text(0.95, 1, reg_func, c="r")
ax.text(0, 0.98, f"PCC: {pearson:.2f}")
ax.set_xlabel("Replicate 1")
ax.set_ylabel("Replicate 2")
ax.set_title(tp_name[tp])
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/Ai/venn/{tp}-AI-scatter.svg",
    ),
    bbox_inches="tight",
)
