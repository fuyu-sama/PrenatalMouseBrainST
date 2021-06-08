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

ncs_full = {
    "E135A": 10,
    "E135B": 12,
    "E155A": 11,
    "E155B": 11,
    "E165A": 11,
    "E165B": 10,
    "E175A1": 10,
    "E175A2": 18,
    "E175B": 12,
    "P0B": 12,
    "P0A1": 14,
    "P0A2": 10
}

# %% read data
cluster_full_df = pd.DataFrame(columns=["sc3_clusters"])
coor_full_df = pd.DataFrame()
for idx in idx_full:
    cluster_path = Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/results/cluster/SCT-SC3/pattern/{idx}-SC3.csv"
    )
    cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)
    cluster_series = cluster_df[f'sc3_{ncs_full[idx]}_clusters']
    cluster_df = pd.DataFrame(
        [f"{idx}_{i}" for i in cluster_series],
        index=cluster_df.index,
    )
    cluster_df.columns = ["sc3_clusters"]
    cluster_full_df = pd.concat([cluster_full_df, cluster_df])

    coor_path = Path.joinpath(
        Path.home(), f"workspace/mouse-brain-full/coor_df/{idx}-coor.csv")
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    coor_full_df = pd.concat([coor_full_df, coor_df])

cluster_full_df.to_csv(
    Path.joinpath(
        Path.home(),
        f"workspace/mouse-brain-full/results/cluster/SCT-SC3/pattern/full-SC3.csv"
    ))
coor_full_df.to_csv(
    Path.joinpath(Path.home(),
                  "workspace/mouse-brain-full/coor_df/full-coor.csv"))
