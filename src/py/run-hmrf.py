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
import session_info
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image

import sakuyahelp

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})
Image.MAX_IMAGE_PIXELS = None
random_state = 42

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
colors = [
    "#FAEBD7", "#00FFFF", "#FFD700", "#0000FF", "#FF8C00", "#EE82EE",
    "#9ACD32", "#5F9EA0", "#7FFF00", "#7FFFD4", "#6495ED", "#008B8B",
    "#B8860B", "#C0C0C0", "#000080", "#D8BFD8", "#00CED1", "#9400D3",
    "#8E804B", "#0089A7", "#CB1B45", "#FFB6C1", "#00FF00", "#800000",
    "#376B6D", "#D8BFD8", "#F5F5F5", "#D2691E"
]

try:
    idx = sys.argv[1]
    scale_method = sys.argv[2]
    cluster_method = sys.argv[3]
except IndexError:
    idx = "E165A"
    scale_method = "combat"
    cluster_method = "sc3"

# %% read data
he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
scale_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv"
)

he_image = Image.open(he_path)
scale_df = pd.read_csv(scale_path, index_col=0, header=0).T
coor_df = pd.read_csv(coor_path, index_col=0, header=0)
cluster_df = pd.read_csv(cluster_path, index_col=0, header=0)
cluster_df = cluster_df.reindex(index=scale_df.index)

# %% run hmrf
hmrf_df = pd.DataFrame(index=cluster_df.index)
cluster_result = sakuyahelp.cluster.HMRF_cluster(
    expression=scale_df,
    coordinate=coor_df,
    init_cluster=[int(i.split("_")[1]) for i in cluster_df[f"{cluster_method}_clusters"]],
    beta=0.5,
)
hmrf_df[f"{cluster_method}-hmrf_clusters"] = cluster_result
