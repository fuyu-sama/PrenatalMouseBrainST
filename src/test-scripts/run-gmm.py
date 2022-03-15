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
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

import SpaGene

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
    scale_method = "combat"

# %%
tp_full = {
    "E135": ["E135A", "E135B"],
    "E155": ["E155A", "E155B"],
    "E165": ["E165A", "E165B"],
    "E175": ["E175A1", "E175A2", "E175B"],
    "P0": ["P0A1", "P0A2"],
}

for n in [2, 3]:
    ss = []
    for tp in tp_full:
        s = []
        for idx in tp_full[tp]:
            global_moran = pd.read_csv(
                Path.joinpath(
                    WORKDIR,
                    f"results/global_moran/{idx}-{scale_method}-8.csv",
                ),
                index_col=0,
                header=0,
            )

            selected_genes, ax, gmm = SpaGene.utils.filter_i(
                global_moran["I_value"].dropna().to_frame(),
                n_components=n,
            )
            s.append(selected_genes)
            [ss.append(i) for i in selected_genes]
            ax.figure.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/5/{idx}/{scale_method}-{n}-pdf.jpg",
                ))
        fig, ax = plt.subplots()
        if len(s) == 2:
            venn2([set(i) for i in s], tp_full[tp], ax=ax)
        elif len(s) == 3:
            venn3([set(i) for i in s], tp_full[tp], ax=ax)
        ax.set_title(f"{tp} n_components = {n}")
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/1/{tp}-{scale_method}-{n}.jpg",
            ))
    with open(
            Path.joinpath(
                WORKDIR,
                f"results/I-gmm/full-{scale_method}-{n}.csv",
            ), "w") as f:
        for i in set(ss):
            f.write(f"{i}\n")
