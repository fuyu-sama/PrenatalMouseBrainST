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
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

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
                 "E175B_11", "P0A1_8", "P0A2_2", "P0B_1"),
)

# %% read data
coor_path = Path.joinpath(WORKDIR, "Data/coor_df/full-coor.csv")
cluster_path = Path.joinpath(
    WORKDIR,
    "results/cluster/SCT-SC3/pattern/full-SC3.csv",
)
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)

# %% draw
for idx in idx_full:
    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    with Image.open(he_path) as he_image:
        for region in regions:
            region_index = [
                i for i in cluster_df.index
                if cluster_df.loc[i, "sc3_clusters"] in
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
                    f"results/cluster/SCT-SC3/region/{idx}_{region}.jpg",
                ))
            plt.close(fig)
