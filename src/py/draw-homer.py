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

import numpy as np
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

alias = {
    "Scl": "Tal1",
    "Oct6": "Pou3f1",
    "Bmyb": "Mybl2",
    "Eklf": "Klf1",
    "Rest-nrsf": "Prickle1",
    "Snail1": "Snai1",
}
gene_family = {
    "Zic": ["Zic1", "Zic4", "Zic5", "Zic2", "Zic3"],
    "Nfy": ["Nfyc", "Nfyb", "Nfya"],
    "Tead": ["Tead4", "Tead2", "Tead1", "Tead3"],
    "Ebf": ["Ebf4", "Ebf3", "Ebf1", "Ebf2"],
    "Rfx": ["Rfx5", "Rfx1", "Rfx7", "Rfx4", "Rfx2", "Rfx3"],
    "Runx": ["Runx1", "Runx2"],
}

try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
except IndexError:
    scale_method = "combat-qq-logcpm"
    cluster_method = "sc3"

# %% read count data
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
count_df = pd.read_csv(count_path, index_col=0, header=0).T
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

# %% read homer
draw_pvalue = pd.DataFrame()
for region in regions:
    known_path = Path.joinpath(
        WORKDIR,
        f"results/motifResults/{scale_method}-{cluster_method}/{region}/knownResults.txt",
    )
    read_df = pd.read_csv(known_path, index_col=0, header=0, sep="\t")
    read_df = read_df[read_df["P-value"] <= 0.01]
    known_df = pd.DataFrame(-read_df["Log P-value"])
    known_df.index = [i.split("/")[0].split("(")[0] for i in read_df.index]
    known_df.columns = [region]
    known_df = known_df.drop_duplicates()
    draw_pvalue = pd.concat([draw_pvalue, known_df], axis=1)

pvalue_min = draw_pvalue.min().min()
pvalue_max = draw_pvalue.max().max()
draw_pvalue = draw_pvalue.fillna(0.0)
draw_pvalue.index = [i.lower().capitalize() for i in draw_pvalue.index]
draw_pvalue.index = [alias[i] if i in alias else i for i in draw_pvalue.index]

not_in = []
for i in draw_pvalue.index:
    if i not in count_df.columns and i not in gene_family:
        not_in.append(i)
draw_pvalue = draw_pvalue.drop(index=not_in)

# %% build draw_count
draw_count = pd.DataFrame(index=draw_pvalue.index, columns=draw_pvalue.columns)
for region in regions:
    in_region = [
        i for i in cluster_df.index
        if cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region]
    ]
    count_sub = count_df.reindex(index=in_region)
    for gene in draw_pvalue.index:
        if gene in gene_family:
            flag = count_sub.reindex(columns=gene_family[gene]).T
            draw_count.loc[gene, region] = flag.mean().mean()
        else:
            draw_count.loc[gene, region] = count_sub[gene].mean()

# %% draw
fig, ax = plt.subplots(figsize=(10, 30))
vmax = draw_count.max().max() * 0.9
vmin = draw_count.min().min()
for i in range(draw_pvalue.shape[1]):
    for j in range(draw_pvalue.shape[0]):
        sc = ax.scatter(
            i,
            j,
            s=draw_pvalue.iloc[j, i] * 50,
            c=draw_count.iloc[j, i],
            cmap="Reds",
            vmax=vmax,
            vmin=vmin,
        )
ax.set_xticks(range(draw_pvalue.shape[1]))
ax.set_yticks(range(draw_pvalue.shape[0]))
ax.set_xticklabels(draw_pvalue.columns, rotation=45, ha="right")
ax.set_yticklabels(draw_pvalue.index)
ax.yaxis.tick_right()
cb = fig.colorbar(sc, location="top")
cb.set_ticks([vmin, vmax])
cb.set_ticklabels(["Low", "High"])
ax_legend = fig.add_axes([0, 0.11, 0.1, 0.2])
ax_legend.scatter(
    [1, 1, 1],
    [0.9, 1, 1.1],
    s=[pvalue_min * 50, -np.log(1e-5) * 50, pvalue_max * 50],
)
ax_legend.set_xticks([])
ax_legend.set_yticks([0.85, 0.9, 1, 1.1, 1.15])
ax_legend.set_yticklabels(["", "1e-2", "1e-5", "1e-7", ""])
ax_legend.set_title("pvalue")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/motifResults/{scale_method}-{cluster_method}/motifResults.jpg",
    ),
    bbox_inches="tight",
)
