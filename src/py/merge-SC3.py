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
from pathlib import Path

import session_info
import pandas as pd

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()

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

scale_method = "combat"

# %% read data
cluster_full_df = pd.DataFrame(columns=["sc3_clusters"])
regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-SC3/regions.json")
with open(regions_path) as f:
    ncs_full = json.load(f)["ncs"]

for idx in idx_full:
    cluster_path = Path.joinpath(
        WORKDIR, f"results/cluster/{scale_method}-SC3/pattern/{idx}-SC3.csv")
    cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)
    cluster_series = cluster_df[f'sc3_{ncs_full[idx]}_clusters']
    cluster_df = pd.DataFrame(
        [f"{idx}_{i}" for i in cluster_series],
        index=cluster_df.index,
    )
    cluster_df.columns = ["sc3_clusters"]
    cluster_full_df = pd.concat([cluster_full_df, cluster_df])

cluster_full_df.to_csv(
    Path.joinpath(
        WORKDIR, f"results/cluster/{scale_method}-SC3/pattern/full-SC3.csv"), )