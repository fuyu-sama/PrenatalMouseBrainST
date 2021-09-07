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

try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
except IndexError:
    scale_method = "combat"
    cluster_method = "sc3"

# %% read data
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]

# %% draw
for idx in idx_full:
    coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    with Image.open(he_path) as he_image:
        for region in regions:
            region_index = [
                i for i in cluster_df.index
                if cluster_df.loc[i, f"{cluster_method}_clusters"] in
                [p for p in regions[region] if idx in p]
            ]
            if len(region_index) == 0:
                continue
            draw_df = coor_df.reindex(index=region_index)
            fig, ax = plt.subplots(figsize=(10, 10))
            ax.imshow(he_image)
            ax.scatter(draw_df['X'], draw_df['Y'], s=16, alpha=0.7)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f'{idx} {region}')
            [axis.set_visible(False) for axis in ax.spines.values()]
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster/{scale_method}-{cluster_method}/region/{idx}_{region}.jpg",
                ))
            plt.close(fig)
