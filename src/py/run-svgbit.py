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
import os
import sys
from pathlib import Path

import pandas as pd
import psutil
from matplotlib import pyplot as plt
from PIL import Image

import svgbit

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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

try:
    idx = sys.argv[1]
except IndexError:
    idx = "E165A"


def printmen():
    print(psutil.Process(os.getpid()).memory_info().rss / 1024**2)


# %% read data
count_path = Path.joinpath(WORKDIR, f"Data/scale_df/raw/{idx}-raw.csv")
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = Image.open(he_path)

write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{idx}")
if not write_dir.exists():
    write_dir.mkdir()

# %%
# for var in [0, 1, 2, 3, 4, 5]:
for var in [3, 4, 5]:
    for exp in [0.95, 0.9, 0.85, 0.8]:
        print(f"var: {var}, exp: {exp}")
        printmen()
        dataset = svgbit.STDataset(count_df, coor_df)
        dataset = svgbit.filters.low_variance_filter(dataset, var)
        dataset = svgbit.filters.high_expression_filter(dataset, exp)
        dataset = svgbit.normalizers.logcpm_normalizer(dataset)
        svgbit.run(dataset, cores=10)
        svgbit.svg_heatmap(
            dataset,
            Path.joinpath(
                WORKDIR,
                f"results/svgbit/{idx}/{idx}-{var}-{exp}.jpg",
            ),
            he_image,
        )
        dataset.svg_cluster.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/svgbit/{idx}/{idx}-{var}-{exp}.csv",
            ))
