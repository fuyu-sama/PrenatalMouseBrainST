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
import session_info
from matplotlib import pyplot as plt
from PIL import Image

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
    # "E175B": "V10M17-085-E175B",
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

idx = sys.argv[1]

# %% read data
count_df = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"Data/scale_df/logcpm-hotspot-6/{idx}-logcpm-hotspot-6.csv",
    ),
    index_col=0,
    header=0,
).T

coor_df = pd.read_csv(
    Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv"),
    index_col=0,
    header=0,
)

he_image = Image.open(Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif"))

# %%
fig, ax = plt.subplots(figsize=(10, 10), dpi=50)
ax.imshow(he_image)
ax.axis("off")
hotspot_list = []
for i, gene in enumerate(["Satb2", "Nr4a2"]):
    draw_counts = count_df[gene].loc[count_df[gene] > 0]
    ax.scatter(
        coor_df["X"].reindex(index=draw_counts.index),
        coor_df["Y"].reindex(index=draw_counts.index),
        c=colors[i + 1],
        alpha=0.7,
        s=16,
        label=f"{gene}",
    )
    hotspot_list.append(draw_counts.index)

inter_spots = set()
for j, k in enumerate(hotspot_list):
    if j == 0:
        inter_spots = set(k)
    else:
        inter_spots = inter_spots & set(k)
ax.scatter(
    coor_df["X"].reindex(index=inter_spots),
    coor_df["Y"].reindex(index=inter_spots),
    c=colors[i + 2],
    alpha=0.7,
    s=16,
    label=f"Both",
)

ax.legend(markerscale=5)
fig.savefig(Path.joinpath(
    WORKDIR,
    f"draw_genes/{idx}.jpg",
))
plt.close(fig)
