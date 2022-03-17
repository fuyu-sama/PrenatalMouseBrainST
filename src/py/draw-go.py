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
    scale_method = "combat-gmm-2"
    cluster_method = "sc3"

# %% read data
regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

pvalue_df = pd.DataFrame()
fraction_df = pd.DataFrame()
for region in regions:
    homer_path = Path.joinpath(
        WORKDIR,
        f"results/motifResults/{scale_method}-{cluster_method}/{region}")
    processes = [
        "biological_process", "molecular_function", "cellular_component"
    ]
    go_df = pd.DataFrame()
    for process in processes:
        read_df = pd.read_csv(
            Path.joinpath(homer_path, f"{process}.txt"),
            index_col=1,
            header=0,
            sep="\t",
        )
        read_df.index = [f"[{process}] {i}" for i in read_df.index]
        read_df = read_df[read_df["logP"] < -3]
        go_df = pd.concat([go_df, read_df.iloc[0:5, ]])
    read_pvalue = pd.DataFrame(-go_df["logP"])
    read_fraction = pd.DataFrame(go_df["Fraction of Targets in Term"])
    read_pvalue.columns = [region]
    read_fraction.columns = [region]
    pvalue_df = pd.concat([pvalue_df, read_pvalue], axis=1)
    fraction_df = pd.concat([fraction_df, read_fraction], axis=1)
pvalue_df = pvalue_df.fillna(0)
fraction_df = fraction_df.fillna(0)
pvalue_df = pvalue_df.sort_index()
fraction_df = fraction_df.sort_index()

# %% draw
fig, ax = plt.subplots(figsize=(15, 45))
vmax = pvalue_df.max().max()
vmin = pvalue_df.min().min()
for i in range(pvalue_df.shape[1]):
    for j in range(pvalue_df.shape[0]):
        sc = ax.scatter(
            i,
            j,
            s=fraction_df.iloc[j, i] * 800,
            c=pvalue_df.iloc[j, i],
            cmap="Reds",
            vmax=vmax,
            vmin=vmin,
        )
ax.set_xticks(range(pvalue_df.shape[1]))
ax.set_yticks(range(pvalue_df.shape[0]))
ax.set_xticklabels(pvalue_df.columns, rotation=45, ha="right")
ax.set_yticklabels(pvalue_df.index)
ax.yaxis.tick_right()
ax_cb = fig.add_axes([0, 0.7, 0.05, 0.1])
cb = fig.colorbar(sc, cax=ax_cb)
cb.set_ticks([vmin, vmax])
cb.set_ticklabels(["Low", "High"])
cb.ax.yaxis.tick_left()
cb.ax.yaxis.set_label_position("left")
cb.ax.set_ylabel("-logP")
ax_legend = fig.add_axes([0, 0.11, 0.1, 0.2])
ax_legend.scatter(
    [1, 1, 1, 1],
    [0.9, 1, 1.1, 1.2],
    s=[
        fraction_df.max().max() * 200,
        fraction_df.max().max() * 400,
        fraction_df.max().max() * 600,
        fraction_df.max().max() * 800
    ],
)
ax_legend.set_xticks([])
ax_legend.set_yticks([0.85, 0.9, 1, 1.1, 1.2, 1.25])
ax_legend.set_yticklabels([
    "", f"{fraction_df.max().max() / 4:.3f}",
    f"{fraction_df.max().max() / 2:.3f}",
    f"{fraction_df.max().max() / 4 * 3:.3f}", f"{fraction_df.max().max():.3f}",
    ""
])
ax_legend.set_title("Gene ratio")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/motifResults/{scale_method}-{cluster_method}/go.jpg"),
    bbox_inches="tight",
)
plt.close(fig)
