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
import seaborn as sns
import matplotlib.pyplot as plt

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
    idx = sys.argv[1]
except IndexError:
    idx = "E165A"

# %%
global_moran_path = Path.joinpath(
    WORKDIR,
    f"results/global_moran/{idx}-logcpm-6.csv",
)
global_moran_df = pd.read_csv(
    global_moran_path,
    index_col=0,
    header=0,
).dropna(axis="index")
global_moran_df = global_moran_df.sort_values(by="I_value", ascending=False)

# %%
draw_list = []
for i in range(0, 5000, 500):
    sub_series = global_moran_df.iloc[i: i + 500, 1]
    sub_series.name = f"{i + 1} - {i + 500}"
    draw_list.append(sub_series)
fig, ax = plt.subplots(figsize=(10, 10))
sns.boxplot(data=draw_list, ax=ax)
ax.set_xticklabels([i.name for i in draw_list], rotation=45)
ax.set_title(f"{idx} Global Moran's I")
ax.set_xlabel("Gene ranking")
ax.set_ylabel("Global Moran's I")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/global_moran/{idx}-logcpm-6.jpg",
    ),
    bbox_inches="tight",
)
