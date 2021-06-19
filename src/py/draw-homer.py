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
regions = dict(
    cortex=("E135A_1", "E135B_3", "E135B_9", "E155A_4", "E155B_5", "E155B_6",
            "E165A_1", "E165B_4", "E175A1_8", "E175A2_5", "E175A2_6",
            "E175B_6", "P0A1_5", "P0A2_1", "P0B_3"),
    thalamus=("E135A_2", "E135B_6", "E135B_7", "E155A_5", "E155A_7", "E155A_8",
              "E155B_7", "E165A_3", "E165A_5", "E165B_3", "E175A1_2",
              "E175A1_7", "E175A2_10", "E175A2_11", "E175A2_12", "E175B_2",
              "E175B_4", "E175B_5", "P0A1_6", "P0A1_7", "P0A1_12", "P0A2_6",
              "P0B_1", "P0B_7", "P0B_8"),
    hypothalamus=("E135A_3", "E135A_6", "E135B_5", "E135B_12", "E155A_6",
                  "E155B_1", "E165A_6", "E165B_1", "E175B_3", "P0A1_1",
                  "P0A2_4", "P0B_12"),
    olfactory=("E135A_5", "E135A_7", "E135B_4", "E155A_11", "E155B_2",
               "E155B_11", "E165A_10", "E165A_11", "E165B_7", "E165B_10",
               "E175A1_5", "E175A2_7", "E175A2_9", "E175B_1", "E175B_9",
               "P0A1_10", "P0A2_7", "P0B_5", "P0B_6"),
    hippocampus=("E155A_2", "E165A_2", "E165B_6", "E175A2_14", "E175A2_15",
                 "E175B_11", "P0A1_8", "P0A2_2", "P0B_4"),
)

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

# %% read count data
count_path = Path.joinpath(
    WORKDIR,
    "Data/scale_df/logcpm/full-logcpm-inter.csv",
)
cluster_path = Path.joinpath(
    WORKDIR,
    "results/cluster/SCT-SC3/pattern/full-SC3.csv",
)
count_df = pd.read_csv(count_path, index_col=0, header=0).T
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)

# %% read homer
draw_pvalue = pd.DataFrame()
for region in regions:
    known_path = Path.joinpath(
        WORKDIR, f"results/motifResults/{region}/knownResults.txt")
    read_df = pd.read_csv(known_path, index_col=0, header=0, sep="\t")
    read_df = read_df[read_df["P-value"] <= 0.01]
    known_df = pd.DataFrame(-read_df["Log P-value"])
    known_df.index = [i.split("/")[0].split("(")[0] for i in read_df.index]
    known_df.columns = [region]
    draw_pvalue = pd.concat([draw_pvalue, known_df], axis=1)
pvalue_min = draw_pvalue.min().min()
pvalue_max = draw_pvalue.max().max()
draw_pvalue = draw_pvalue.fillna(0.5)
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
        if cluster_df.loc[i, "sc3_clusters"] in regions[region]
    ]
    count_sub = count_df.reindex(index=in_region)
    for gene in draw_pvalue.index:
        if gene in gene_family:
            flag = count_sub.reindex(columns=gene_family[gene]).T
            draw_count.loc[gene, region] = flag.mean().mean()
        else:
            draw_count.loc[gene, region] = count_sub[gene].mean()

# %% draw
fig, ax = plt.subplots(figsize=(10, 25))
vmax = draw_count.max().max()
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
ax.set_xticklabels(draw_pvalue.columns)
ax.set_yticklabels(draw_pvalue.index)
fig.colorbar(sc, location="top")
ax_legend = fig.add_axes([1, 0.11, 0.1, 0.2])
ax_legend.scatter(
    [1, 1, 1],
    [0.9, 1, 1.1],
    s=[0.5 * 50, pvalue_min * 50, pvalue_max * 50],
)
ax_legend.set_xticks([])
ax_legend.set_yticks([0.85, 0.9, 1, 1.1, 1.15])
ax_legend.set_yticklabels(["", "1", "1e-2", "1e-7", ""])
ax_legend.yaxis.tick_right()
ax_legend.set_title("pvalue")
fig.savefig(
    Path.joinpath(WORKDIR, "results/motifResults/motifResults.jpg"),
    bbox_inches="tight",
)
