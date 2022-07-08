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
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

# %% read data
idx = "E165A"

count_df = pd.read_csv(
    Path.joinpath(WORKDIR, f"Data/scale_df/logcpm/{idx}-logcpm.csv"),
    index_col=0,
    header=0,
).T

hotspot_df = pd.read_csv(
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

ai_df = pd.read_csv(
    Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
    index_col=0,
    header=0,
).sort_values(by="Ai")

he_image = Image.open(Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif"))

# %% draw
for gene in ai_df[ai_df["Ai"] <= 0].index:
    if gene[:2] == "Gm":
        continue
    if gene[-3:] == "Rik":
        continue
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(26, 10))
    ax1.imshow(he_image)
    ax1.axis("off")
    ax1.set_title(f"{gene} logcpm")
    sc = ax1.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=count_df[gene],
        cmap='autumn_r',
        alpha=0.7,
        s=16,
    )
    cb = fig.colorbar(sc, ax=ax1)

    ax2.imshow(he_image)
    ax2.axis("off")
    ax2.set_title(f"{gene} hotspot, AI = {ai_df.loc[gene, 'Ai']:.3f}")
    sc = ax2.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=hotspot_df[gene],
        cmap='autumn_r',
        alpha=0.7,
        s=16,
        vmin=0,
        vmax=1,
    )
    cb = fig.colorbar(sc, ax=ax2)
    fig.savefig(Path.joinpath(
        WORKDIR,
        f"draw_genes/1/{idx}-{gene}.jpg",
    ))
    plt.close(fig)
