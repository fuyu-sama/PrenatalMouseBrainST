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
    scale_method = "combat-Ai-2000_2500"

# %% read data
with open(
        Path.joinpath(
            WORKDIR,
            f"results/cluster/{scale_method}-sc3/color_swift.json",
        )) as f:
    color_swift = json.load(f)
for idx, colors in color_swift.items():
    ncs = len(colors)
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
    for i in colors:
        draw_series.replace(int(i), colors[i], inplace=True)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=50)
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