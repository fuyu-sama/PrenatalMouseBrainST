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

import svgbit

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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

try:
    idx = sys.argv[1]
except IndexError:
    idx = "E165A"

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

# %%
for q in [0.95, 0.99]:
    dataset = svgbit.STDataset(count_df, coor_df)
    dataset = svgbit.filters.low_variance_filter(dataset, 0)
    dataset_new = svgbit.filters.quantile_filter(dataset, q)
    dataset_new = svgbit.normalizers.logcpm_normalizer(dataset_new)
    svgbit.run(dataset_new, cores=5, n_svg_clusters=9)
    dataset_new.AI.sort_values(ascending=False).to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/svgbit/{idx}-{q}-AI.csv",
        ))
    dataset_new.svg_cluster.to_csv(
        Path.joinpath(
            WORKDIR,
            f"results/svgbit/{idx}-{q}-cluster.csv",
        ))
    svgbit.svg_heatmap(
        dataset_new,
        Path.joinpath(WORKDIR, f"results/svgbit/{idx}-{q}.jpg"),
        he_path,
    )
    with open(Path.joinpath(WORKDIR, f"results/svgbit/{idx}-{q}-drop.csv"),
              "w") as f:
        for i in dataset.genes:
            if i not in dataset_new.genes:
                f.write(f"{i}\n")
