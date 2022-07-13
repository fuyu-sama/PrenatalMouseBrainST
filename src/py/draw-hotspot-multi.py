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
from itertools import combinations
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

# %% read data
count_dict = {}
coor_dict = {}
he_dict = {}
for idx in idx_full:
    count_dict[idx] = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"Data/scale_df/logcpm-hotspot-6/{idx}-logcpm-hotspot-6.csv",
        ),
        index_col=0,
        header=0,
    ).T

    coor_dict[idx] = pd.read_csv(
        Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv"),
        index_col=0,
        header=0,
    )

    he_dict[idx] = Image.open(
        Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif"))

# %%
colors = [
    "tab:red",
    "tab:green",
    "tab:blue",
    "tab:orange",
    "tab:pink",
    "tab:cyan",
    "k",
]
genes = [
    "Satb2", "Tbr1", "Eomes"
]

fig, axes = plt.subplots(2, 5, figsize=(15 * 5, 20), dpi=500)
for idx, ax in zip(idx_full.keys(), axes.flatten()):
    ax.imshow(he_dict[idx])
    ax.axis("off")
    for i, gene in enumerate(genes):
        draw_counts = count_dict[idx][gene].loc[count_dict[idx][gene] == 1]
        ax.scatter(
            coor_dict[idx]["X"].reindex(index=draw_counts.index),
            coor_dict[idx]["Y"].reindex(index=draw_counts.index),
            c=colors[i],
            alpha=0.7,
            s=16,
            label=f"{gene}",
        )

    n_combinations = 2
    while n_combinations <= len(genes):
        for c in combinations(genes, n_combinations):
            i += 1
            draw_counts = pd.Series(index=coor_dict[idx].index, dtype=int)
            draw_counts = draw_counts.fillna(0)
            for gene in c:
                draw_counts += count_dict[idx][gene]
            draw_counts = draw_counts[draw_counts == n_combinations]
            ax.scatter(
                coor_dict[idx]["X"].reindex(index=draw_counts.index),
                coor_dict[idx]["Y"].reindex(index=draw_counts.index),
                c=colors[i],
                s=16,
                label=f"{' & '.join(c)}",
            )
        n_combinations += 1
    ax.legend(
        markerscale=2,
        bbox_to_anchor=(1.75, 1),
        loc='upper right',
    )
    ax.set_title(idx)

fig.savefig(Path.joinpath(
    WORKDIR,
    f"draw_genes/{'-'.join(genes)}-high.jpg",
), )
plt.close(fig)

fig, axes = plt.subplots(2, 5, figsize=(15 * 5, 20), dpi=200)
for idx, ax in zip(idx_full.keys(), axes.flatten()):
    ax.imshow(he_dict[idx])
    ax.axis("off")
    for i, gene in enumerate(genes):
        draw_counts = count_dict[idx][gene].loc[count_dict[idx][gene] == 1]
        ax.scatter(
            coor_dict[idx]["X"].reindex(index=draw_counts.index),
            coor_dict[idx]["Y"].reindex(index=draw_counts.index),
            c=colors[i],
            alpha=0.7,
            s=16,
            label=f"{gene}",
        )

    n_combinations = 2
    while n_combinations <= len(genes):
        for c in combinations(genes, n_combinations):
            i += 1
            draw_counts = pd.Series(index=coor_dict[idx].index, dtype=int)
            draw_counts = draw_counts.fillna(0)
            for gene in c:
                draw_counts += count_dict[idx][gene]
            draw_counts = draw_counts[draw_counts == n_combinations]
            ax.scatter(
                coor_dict[idx]["X"].reindex(index=draw_counts.index),
                coor_dict[idx]["Y"].reindex(index=draw_counts.index),
                c=colors[i],
                s=16,
                label=f"{' & '.join(c)}",
            )
        n_combinations += 1
    ax.legend(
        markerscale=2,
        bbox_to_anchor=(1.75, 1),
        loc='upper right',
    )
    ax.set_title(idx)

fig.savefig(Path.joinpath(
    WORKDIR,
    f"draw_genes/{'-'.join(genes)}-low.jpg",
), )
plt.close(fig)
