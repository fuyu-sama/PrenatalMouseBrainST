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

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3

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
    "E175B": "E175",
    # "P0B": "P0",
    "P0A1": "P0",
    "P0A2": "P0",
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
    tp = "E165"

# %%
moran_list = []
for idx in tp_full[tp]:
    moran_list.append(
        pd.read_csv(
            Path.joinpath(WORKDIR, f"results/global_moran/{idx}-logcpm-6.csv"),
            index_col=0,
            header=0,
        ).sort_values(by="I_value", ascending=False))

# %%
len_inter = []
venn_function = venn2 if len(moran_list) == 2 else venn3
for i in range(0, 5000, 500):
    fig, ax = plt.subplots(figsize=(10, 10))
    set_list = [set(j.iloc[i:i + 500, 0]) for j in moran_list]
    venn_function(set_list, tp_full[tp], ax=ax)
    ax.set_title(f"{tp} {i} - {i + 500}")
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/global_moran/venn/{tp}-{i}_{i + 500}.jpg",
        ))
    plt.close(fig)

    for j, k in enumerate(set_list):
        if j == 0:
            inters = k
        else:
            inters = inters & k
    len_inter.append(len(inters))

fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(range(len(len_inter)), len_inter, "D-")
ax.set_xticks(range(10))
ax.set_xticklabels(
    [f"{i} - {i + 500}" for i in range(0, 5000, 500)],
    rotation=45,
)
ax.set_title(f"{tp} intersection genes")
ax.set_xlabel("Gene ranking")
ax.set_ylabel("Intersection genes")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/global_moran/venn/{tp}-plot.jpg",
    ),
    bbox_inches="tight",
)
