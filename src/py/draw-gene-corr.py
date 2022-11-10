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

# %%
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

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
colors = [
    "#FAEBD7", "#00FFFF", "#FFD700", "#0000FF", "#FF8C00", "#EE82EE",
    "#9ACD32", "#5F9EA0", "#7FFF00", "#7FFFD4", "#6495ED", "#008B8B",
    "#B8860B", "#C0C0C0", "#000080", "#D8BFD8", "#00CED1", "#9400D3",
    "#8E804B", "#0089A7", "#CB1B45", "#FFB6C1", "#00FF00", "#800000",
    "#376B6D", "#D8BFD8", "#F5F5F5", "#D2691E"
]

try:
    idx = sys.argv[1]
    scale_method = sys.argv[2]
except IndexError:
    idx = "E175A1"
    scale_method = "logcpm"

# %%
count_df = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
    ),
    index_col=0,
    header=0,
).T

with open(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster/logcpm-hotspot-6-logcpm-0.99-Ai-0_500/regions.json"
        )) as f:
    regions = json.load(f)
ncs = len(regions[idx])

cluster_series = pd.read_csv(
    Path.joinpath(
        WORKDIR,
        f"results/gene-cluster/logcpm-hotspot-6-logcpm-0.99-Ai-0_500",
        f"{idx}-{ncs}/{idx}-spots.csv",
    ),
    index_col=0,
    header=0,
)[f"type_1"]
cortex_cluster = [int(i) for i in regions[idx] if "cortex" in regions[idx][i]]
hippo_cluster = [
    int(i) for i in regions[idx] if "hippocampus" in regions[idx][i]
]

# %%
gene_1 = ["Satb2", "Mef2c"]
gene_2 = ["Nr4a2", "Mef2c", "Bcl11b", "Dkk3", "Zbtb20"]
gene_1 = ["Satb2"]
gene_2 = ["Fezf2"]
for q in gene_1:
    for p in gene_2:
        fig, axes = plt.subplots(1, 3, figsize=(30, 10))
        draw_spots = cluster_series.index
        axes[0].scatter(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )
        corr = np.corrcoef(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )[1, 0]
        axes[0].set_title(f"All spots\npearson = {corr:.4f}")
        draw_spots = [
            j for i in cortex_cluster
            for j in cluster_series[cluster_series == i].index
        ]
        axes[1].scatter(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )
        corr = np.corrcoef(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )[1, 0]
        axes[1].set_title(f"cortex spots\npearson = {corr:.4f}")
        draw_spots = [
            j for i in hippo_cluster + cortex_cluster
            for j in cluster_series[cluster_series == i].index
        ]
        corr = np.corrcoef(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )[1, 0]
        axes[2].scatter(
            count_df[q].reindex(index=draw_spots),
            count_df[p].reindex(index=draw_spots),
        )
        axes[2].set_title(f"cortex + hippocampus\npearson = {corr:.4f}")
        [ax.set_xlabel(f"{q} {scale_method}") for ax in axes]
        [ax.set_ylabel(f"{p} {scale_method}") for ax in axes]
        fig.savefig(Path.joinpath(WORKDIR, f"draw_genes/1/{q}-{p}.jpg"))
