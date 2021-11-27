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

# %% envoronment config
import sys
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

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
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
c = {
    "E135": "#00FFFF",
    "E155": "#FFD700",
    "E165": "#EE82EE",
    "E175": "#9ACD32",
    "P0": "#9400D3",
}
try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
except IndexError:
    scale_method = "combat-qq-logcpm"
    cluster_method = "sc3"

# %% read scenic
auc_path = Path.joinpath(
    WORKDIR,
    f"results/SCENIC/{scale_method}-{cluster_method}/auc.csv",
)
auc_df = pd.read_csv(
    auc_path,
    index_col=0,
    header=0,
)

# %% PCA
region = "hippocampus"
fit_df = auc_df.reindex(index=[i for i in auc_df.index if region in i])
result_pca = PCA(n_components=2).fit_transform(fit_df)

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(result_pca[:, 0], result_pca[:, 1])
for i, j in enumerate(fit_df.index):
    ax.annotate(j.split("_"), (result_pca[i, 0], result_pca[i, 1]))
ax.set_title(region)
fig.savefig(f"results/{region}.jpg")
