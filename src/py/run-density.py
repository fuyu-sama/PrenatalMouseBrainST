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
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from libpysal.weights import KNN

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
    "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
tp_full = {
    "E135": ["E135A", "E135B"],
    "E155": ["E155A", "E155B"],
    "E165": ["E165A", "E165B"],
    "E175": ["E175A1", "E175A2"],
    # "E175": ["E175A1", "E175A2", "E175B"],
    "P0": ["P0A1", "P0A2"],
}
idx_tp = {
    "E135A": "E135",
    "E135B": "E135",
    "E155A": "E155",
    "E155B": "E155",
    "E165A": "E165",
    "E165B": "E165",
    "E175A1": "E175",
    "E175A2": "E175",
    "E175B": "E175",
    # "P0B": "P0",
    "P0A1": "P0",
    "P0A2": "P0",
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
except IndexError:
    idx = "E165A"


def hotspot_density_single(gene, hotspot_df, knn_df, weight_df=None):
    gene_hotspot_df = hotspot_df[hotspot_df[gene] != 0]
    gene_hotspot_df = pd.DataFrame(gene_hotspot_df[gene], columns=[gene])
    n_hotspots = len(gene_hotspot_df)
    hotspots = gene_hotspot_df.index.tolist()

    hs = {}
    for i in hotspots:
        hi_wnn_df = pd.DataFrame(knn_df[i], columns=[i])
        hi_wnn_df = pd.DataFrame(hi_wnn_df[hi_wnn_df[i] == 1])
        knn_coors = hi_wnn_df.index.tolist()
        inter_spots = set(hotspots).intersection(set(knn_coors))
        if weight_df is None:
            hs[i] = [len(inter_spots) / knn_df[i].sum()]
        else:
            hs[i] = [(weight_df[gene][inter_spots] /
                      weight_df[gene][knn_coors].sum()).sum()]

    dens = np.array([i[0] for i in hs.values()])
    dens_sum = dens.sum()
    dens_mean = dens_sum / n_hotspots
    dens_df = pd.DataFrame([dens_mean], columns=['density'])
    dens_df.index = [gene]
    weighted_hotspot = pd.DataFrame.from_dict(hs).T
    weighted_hotspot.columns = [gene]
    return {"density": dens_df, "weighted_hotspot": weighted_hotspot}


def hotspot_density(gene_list, hotspot_df, knn_df, weight_df=None, cores=4):
    partial_func = partial(
        hotspot_density_single,
        hotspot_df=hotspot_df,
        knn_df=knn_df,
        weight_df=weight_df,
    )
    pool = Pool(processes=cores)
    result_lists = pool.map(partial_func, gene_list)
    pool.close()
    pool.join()
    density_df = pd.concat([i["density"] for i in result_lists], axis=0)
    weighted_hotspot_df = pd.concat(
        [i["weighted_hotspot"] for i in result_lists],
        axis=1,
    ).fillna(0)

    return {"density": density_df, "weighted_hotspot": weighted_hotspot_df}


def acquire_knn_df(coor_df, n):
    points = np.array(coor_df[['X', 'Y']])
    wnn = KNN.from_array(points, n)
    wnn_df = pd.DataFrame(*wnn.full()).astype(int)
    wnn_df.columns = coor_df.index.tolist()
    wnn_df.index = coor_df.index.tolist()
    return wnn_df.T


# %% read data
hotspot_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/logcpm-hotspot-6/{idx}-logcpm-hotspot-6.csv",
)
hotspot_df = pd.read_csv(
    hotspot_path,
    index_col=0,
    header=0,
).T

weight_path = Path.joinpath(
    WORKDIR,
    f"results/local_moran/{idx}-logcpm-6-P.csv",
)
weight_df = pd.read_csv(
    weight_path,
    index_col=0,
    header=0,
).T

global_moran_path = Path.joinpath(
    WORKDIR,
    f"results/global_moran/{idx}-logcpm-6.csv",
)
global_moran_df = pd.read_csv(
    global_moran_path,
    index_col=0,
    header=0,
).dropna(axis="index")

coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

hotspot_sum = hotspot_df.sum()
remain_genes = hotspot_sum[hotspot_sum != 0].index
hotspot_df = hotspot_df.reindex(columns=remain_genes)
weight_df = weight_df.reindex(columns=remain_genes)

# %%
results = hotspot_density(
    gene_list=hotspot_df.columns,
    hotspot_df=hotspot_df,
    knn_df=acquire_knn_df(coor_df, 6),
    weight_df=-np.log(weight_df),
    cores=8,
)
density_df = results["density"]
results["weighted_hotspot"].T.to_csv(
    Path.joinpath(
        WORKDIR,
        f"Data/scale_df/logcpm-hotspot-6-weighted/{idx}-logcpm-hotspot-6-weighted.csv",
    ))

# %%
fig, ax = plt.subplots(figsize=(10, 10))
ax.hist(density_df, bins=100)
ax.set_title(f"{idx} -logP weighted density")
fig.savefig(Path.joinpath(WORKDIR, f"results/density/{idx}-dist.jpg"))
plt.close()

# %%
draw_df = pd.concat(
    [density_df, global_moran_df["I_value"].to_frame()],
    axis="columns",
).dropna(axis="index")
draw_df.to_csv(Path.joinpath(WORKDIR, f"results/density/{idx}-density.csv"))
r = np.corrcoef(draw_df["density"], draw_df["I_value"])[0, 1]
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(draw_df["density"], draw_df["I_value"], s=8)
ax.set_xlabel("density")
ax.set_ylabel("I_value")
ax.set_title(f"{idx} pearson correlation: {r:.4f}")
fig.savefig(Path.joinpath(WORKDIR, f"results/density/{idx}-corr.jpg"))
plt.close()
