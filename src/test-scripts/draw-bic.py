#! /usr/bin/env Rscript

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
from matplotlib import pyplot as plt

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")

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
    idx = sys.argv[3]

except IndexError:
    scale_method = "combat-zjq"
    cluster_method = "sc3"
    idx = "E165A"

# %% read data
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}-moran/{idx}-{scale_method}-moran.csv",
)
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

drop_genes = []
for i in count_df:
    if count_df[i].var() < 1e-6:
        drop_genes.append(i)
count_df.drop(columns=drop_genes, inplace=True)

cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
cluster_df = cluster_df.reindex(
    index=[i for i in cluster_df.index if idx in i])

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

regions_label = {j: i for i, j in enumerate(regions)}
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"{cluster_method}_clusters"]))
    if i not in in_regions
]

# %%
for method in [
        "BCBimax", "BCCC", "BCPlaid", "BCQuest", "BCSpectral", "BCXmotifs"
]:
    call = ""
    s1 = []
    s2 = []
    with open(Path.joinpath(WORKDIR, f"results/4/{method}.csv")) as f:
        flag = 0
        for line in f:
            line = line.strip()
            if flag == 0:
                call = line.replace(" ", "").replace("language", "")
                flag = 1
            elif flag == 1:
                flag = 2
            elif flag == 2:
                [s1.append(i) for i in line.split(",") if i not in s1]
                flag = 3
            elif flag == 3:
                [s2.append(i) for i in line.split(",") if i not in s2]
                flag = 1

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(count_df.reindex(index=s2, columns=s1), cmap="Reds", aspect="auto")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("Gene")
    ax.set_ylabel("Spot")
    # ax.set_xticks(range(len(s1)))
    # ax.set_yticks(range(len(s2)))
    # ax.set_xticklabels(s1, rotation=90)
    # ax.set_yticklabels(s2)
    ax.set_title(
        f"{call[call.find('method'):-1] }\n{len(s1)} genes, {len(s2)} spots")
    fig.savefig(Path.joinpath(WORKDIR, f"results/4/{method}.jpg"))
