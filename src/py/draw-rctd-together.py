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
import json
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
    region = sys.argv[2]
    scale_method = sys.argv[3]
    cluster_method = sys.argv[4]
except IndexError:
    idx = "E175A1"
    region = "cortex"
    scale_method = "combat-qq-logcpm"
    cluster_method = "sc3"

# %% read data
he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = Image.open(he_path)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv"
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]
in_region = [
    i for i in cluster_df.index if (idx in i and cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region])
]

results_path = Path.joinpath(
    WORKDIR, f"results/RCTD/{region}/{idx}/results/results_df.csv")
results_df = pd.read_csv(results_path, index_col=0, header=0)
weights_path = Path.joinpath(
    WORKDIR, f"results/RCTD/{region}/{idx}/results/weights.csv")
weights_df = pd.read_csv(weights_path, index_col=0, header=0)

coor_df = coor_df.reindex(index=in_region)
results_df = results_df.reindex(index=in_region)
weights_df = weights_df.reindex(index=in_region)

# %% draw weight
wanted_types = {
    "Excitatory.neurons": "#D8BFD8",
    "Inhibitory.neurons": "#7FFFD4",
    "Radial.glia.or.Progenitors": "#B8860B",
    "Vascular.endothelial.cells": "#CB1B45",
    "Tanycytes.α": "#6495ED",
    "Tanycytes.β": "#6495ED",
} if region == "hypothalamus" else {
    "Immature.neurons": "#D8BFD8",
    "Migrating.neurons": "#7FFFD4",
    "Interneurons": "#B8860B",
    "CThPN": "#CB1B45",
    "DL.CPN": "#6495ED",
    "Layer.4": "#00FF00",
    "SCPN": "#FFD700",
}
weights_df = weights_df.reindex(columns=wanted_types.keys())
weights_df = weights_df.dropna(axis=1)
fig, ax = plt.subplots(figsize=(15, 10))
ax.imshow(he_image)
ax.axis("off")
ax.set_title(idx)
for i in weights_df.columns:
    draw_spots = []
    for j in weights_df.index:
        if weights_df.loc[j, i] == weights_df.T.max()[j]:
            draw_spots.append(j)
    sc = ax.scatter(
        coor_df["X"].reindex(index=draw_spots),
        coor_df["Y"].reindex(index=draw_spots),
        s=16,
        c=wanted_types[i],
        label=i,
    )
fig.legend(fontsize=16, markerscale=3)
fig.savefig(
    Path.joinpath(WORKDIR, f"results/RCTD/{region}/{idx}/types.jpg"),
    bbox_inches="tight",
)
