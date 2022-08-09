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

write_dir = Path.joinpath(WORKDIR, f"results/svgbit/{idx}")
if not write_dir.exists():
    write_dir.mkdir()

# %%
var = 0
for exp in [0.95, 0.9, 0.85, 0.8, 0.7]:
    for q in [0.95, 0.9, 0.85, 0.8]:
        dataset = svgbit.STDataset(count_df, coor_df)
        dataset = svgbit.filters.low_variance_filter(dataset, var)
        dataset = svgbit.filters.high_expression_filter(dataset, exp)
        dataset = svgbit.filters.quantile_filter(dataset, q)
        dataset = svgbit.normalizers.logcpm_normalizer(dataset)
        print(f"{idx}, exp: {exp}, quantile: {q}, genes: {dataset.n_genes}")
        svgbit.run(dataset, cores=20, n_svg_clusters=9)
        svgbit.svg_heatmap(
            dataset,
            Path.joinpath(
                WORKDIR,
                f"results/svgbit/{idx}/{idx}-{exp}-{q}.jpg",
            ),
            he_path,
        )
        dataset.svg_cluster.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/svgbit/{idx}/{idx}-{exp}-{q}.csv",
            ))
