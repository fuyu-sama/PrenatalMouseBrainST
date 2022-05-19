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
import os
import sys
from pathlib import Path

import pandas as pd
from anytree import Node, RenderTree
from anytree.exporter import DotExporter
from matplotlib import pyplot as plt

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
    "E155": ["E155B", "E155A"],
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
    "P0B": "P0",
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
    scale_method = sys.argv[1]
except IndexError:
    scale_method = "logcpm"

# %% read data
ncs_path = Path.joinpath(WORKDIR, f"results/gene-cluster/gene-cluster.json")
with open(ncs_path) as f:
    ncs_dict = json.load(f)

cluster_dict = {}
for tp in tp_full:
    idx = tp_full[tp][0]
    cluster_dict[tp] = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"results/gene-cluster/logcpm",
            f"{idx}-{ncs_dict[idx]}/{idx}-genes.csv",
        ),
        index_col=0,
        header=0,
    )

# %%
node_list = []
for i, tp in enumerate(cluster_dict):
    idx = tp_full[tp][0]
    begin_tp = 0
    nodes = []
    cluster_df = cluster_dict[tp]
    if i < begin_tp:
        node_list.append(nodes)
        continue
    elif i == begin_tp:
        for j in set(cluster_df["0"]):
            name = f"{idx}-{j}"
            nodes.append(Node(name, tp=tp, idx=idx, cluster=j))
        node_list.append(nodes)
        continue

    for j in set(cluster_df["0"]):
        cluster_genes = set(cluster_df["0"][cluster_df["0"] == j].index)
        max_coverage, max_node = 0, None
        for k in node_list[i - 1]:
            last_df = cluster_dict[k.tp]
            last_genes = set(last_df["0"][last_df["0"] == k.cluster].index)
            covered_genes = list(last_genes & cluster_genes)
            coverage = len(covered_genes) / len(cluster_genes)
            if coverage > max_coverage:
                max_coverage, max_node = coverage, k
        name = f"{idx}-{j}"
        nodes.append(Node(name, parent=max_node, idx=idx, tp=tp, cluster=j))
    node_list.append(nodes)

# %%
for i, j in enumerate(node_list):
    savedir = Path.joinpath(WORKDIR, f"results/trajectory/{i}")
    if not os.path.exists(savedir):
        os.mkdir(savedir)
    for k, l in enumerate(j):
        DotExporter(l).to_picture(
            Path.joinpath(WORKDIR, savedir, f"{k + 1}.jpg"))
