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
    "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
try:
    scale_method = sys.argv[1]
    n = int(sys.argv[2])
except IndexError:
    scale_method = "combat"
    n = 1000

# %% read data
full_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
full_df = pd.read_csv(full_path, index_col=0, header=0)

# %%
gene_list = []
for idx in idx_full:
    moran_df = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/global_moran/{idx}-logcpm-8.csv"),
        header=0,
        index_col=0,
    ).sort_values(by="I_value", ascending=False)
    [gene_list.append(i) for i in moran_df.index[:n]]
gene_list = set(gene_list)

# %% subset and save
full_df = full_df.reindex(index=set(gene_list))
full_df.to_csv(
    Path.joinpath(
        WORKDIR,
        f"Data/scale_df/{scale_method}-logcpm-{n}",
        f"full-{scale_method}-logcpm-{n}.csv",
    ))
for idx in idx_full:
    sub_df = full_df.reindex(columns=[i for i in full_df.columns if idx in i])
    sub_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"Data/scale_df/{scale_method}-logcpm-{n}",
            f"{idx}-{scale_method}-logcpm-{n}.csv",
        ))
