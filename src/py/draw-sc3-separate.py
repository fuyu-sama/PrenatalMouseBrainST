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

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from PIL import Image

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})
Image.MAX_IMAGE_PIXELS = None

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

try:
    scale_method = sys.argv[1]
except IndexError:
    scale_method = "combat-Ai-2000_2500"

# %% read data
with open(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/regions.json",
        )) as f:
    cluster_regions = json.load(f)
with open(Path.joinpath(WORKDIR, f"results/cluster/colors.json")) as f:
    colors = json.load(f)
for idx, regions in cluster_regions.items():
    ncs = len(regions)
    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    sc3_path = Path.joinpath(
        WORKDIR,
        f"results/cluster/{scale_method}-sc3/pattern/{idx}-sc3.csv",
    )
    coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")

    he_image = Image.open(he_path)
    cluster_df = pd.read_csv(sc3_path, index_col=0, header=0)
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    assert all(cluster_df.index == coor_df.index)

    draw_series = cluster_df[f"sc3_{ncs}_clusters"]
    for i in regions:
        if regions[i] in colors:
            flag = colors[regions[i]]
        else:
            flag = regions[i]
        draw_series.replace(int(i), flag, inplace=True)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
    ax.axis("off")
    ax.imshow(he_image)
    ax.set_title(f"{idx} ncs {ncs}")
    ax.scatter(
        coor_df["X"],
        coor_df["Y"],
        s=16,
        c=cluster_df[f"sc3_{ncs}_clusters"],
        alpha=0.7,
    )
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/{idx}-final.jpg",
        ))
    plt.close(fig)
