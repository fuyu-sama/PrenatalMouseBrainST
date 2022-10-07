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
import random
import sys
from pathlib import Path

import anndata
import cv2
import numpy as np
import pandas as pd
import scanpy as sc
import SpaGCN as spg
import torch
from matplotlib import pyplot as plt

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
    idx = "E135B"

# %%
count_path = Path.joinpath(WORKDIR, f"Data/scale_df/raw/{idx}-raw.csv")
count_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T
count_df = count_df.T[count_df.sum(0) >= 3].T

coor_path = Path.joinpath(
    WORKDIR,
    f"spaceranger/{idx}/outs/spatial/tissue_positions_list.csv",
)
coor_df = pd.read_csv(coor_path, index_col=0, header=None)
coor_df.index = [f"{idx}_{i}" for i in coor_df.index]
coor_df = coor_df.reindex(index=count_df.index)

he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = cv2.imread(str(he_path))

adata = anndata.AnnData(count_df)
adata.obs["x_pixel"] = coor_df[4]
adata.obs["y_pixel"] = coor_df[5]

write_dir = Path.joinpath(WORKDIR, f"results/spagcn/{idx}")
if not write_dir.exists():
    write_dir.mkdir()

# %%
b = 49  # determines the area of each spot when extracting color intensity
b_half = round(b / 2)
x_pixel = adata.obs["x_pixel"].tolist()
y_pixel = adata.obs["y_pixel"].tolist()
x_array = coor_df[2]
y_array = coor_df[3]
img_new = he_image.copy()
for i in range(len(x_pixel)):
    x = x_pixel[i]
    y = y_pixel[i]
    img_new[int(x - b_half):int(x + b_half),
            int(y - b_half):int(y + b_half), :] = 0
img_new[0:100, 0:200, :] = 0
cv2.imwrite(str(Path.joinpath(write_dir, f"{idx}_map.jpg")), img_new)

# %%
adj = spg.calculate_adj_matrix(
    x=x_pixel,
    y=y_pixel,
    x_pixel=x_pixel,
    y_pixel=y_pixel,
    image=he_image,
    beta=b,
    alpha=1,  # determines the weight given to histology
    histology=True,
)

# %%
spg.prefilter_genes(adata, min_cells=3)
spg.prefilter_specialgenes(adata)
sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)

# %%
p = 0.5
l_value = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

n_clusters = 12
r_seed = t_seed = n_seed = 100
res = spg.search_res(
    adata,
    adj,
    l_value,
    n_clusters,
    start=0.7,
    step=0.1,
    tol=5e-3,
    lr=0.05,
    max_epochs=20,
    r_seed=r_seed,
    t_seed=t_seed,
    n_seed=n_seed,
)

# %%
clf = spg.SpaGCN()
clf.set_l(l_value)

random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)

clf.train(
    adata,
    adj,
    init_spa=True,
    init="louvain",
    res=res,
    tol=5e-3,
    lr=0.05,
    max_epochs=200,
)
y_pred, prob = clf.predict()
adata.obs["pred"] = y_pred
adata.obs["pred"] = adata.obs["pred"].astype('category')

# %%
# Do cluster refinement(optional)
# shape="hexagon" for Visium data, "square" for ST data.
adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
refined_pred = spg.refine(
    sample_id=adata.obs.index.tolist(),
    pred=adata.obs["pred"].tolist(),
    dis=adj_2d,
    shape="hexagon",
)
adata.obs["refined_pred"] = refined_pred
adata.obs["refined_pred"] = adata.obs["refined_pred"].astype('category')

# %%
fig, axes = plt.subplots(1, 2, figsize=(20, 10))
for ax in axes:
    ax.imshow(he_image)
    ax.axis("off")
    sc = ax.scatter(y_pixel, x_pixel, s=16, c=y_pred, cmap="tab20")
    leg = ax.legend(
        *sc.legend_elements(),
        ncol=2,
        title="Cluster",
    )
    ax.add_artist(leg)
axes[0].set_title("SpaGCN pred")
axes[1].set_title("SpaGCN refined pred")
fig.savefig(Path.joinpath(write_dir, f"{idx}_pred.jpg"))

# %%
de_genes_info_full = pd.DataFrame()
for target in range(n_clusters):
    adj_2d = spg.calculate_adj_matrix(x=x_array, y=y_array, histology=False)
    start = np.quantile(adj_2d[adj_2d != 0], q=0.001)
    end = np.quantile(adj_2d[adj_2d != 0], q=0.1)
    r = spg.search_radius(
        target_cluster=target,
        cell_id=adata.obs.index.tolist(),
        x=x_array,
        y=y_array,
        pred=adata.obs["pred"].tolist(),
        start=start,
        end=end,
        num_min=10,
        num_max=14,
        max_run=100,
    )
    # Detect neighboring domains
    nbr_domians = spg.find_neighbor_clusters(
        target_cluster=target,
        cell_id=adata.obs.index.tolist(),
        x=x_array,
        y=y_array,
        pred=adata.obs["pred"].tolist(),
        radius=r,
        ratio=1 / 2,
    )

    if nbr_domians is None:
        continue
    nbr_domians = nbr_domians[0:3]
    de_genes_info = spg.rank_genes_groups(
        input_adata=adata,
        target_cluster=target,
        nbr_list=nbr_domians,
        label_col="pred",
        adj_nbr=True,
        log=True,
    )
    de_genes_info["cluster"] = [target] * len(de_genes_info)
    de_genes_info_full = pd.concat([de_genes_info_full, de_genes_info])
de_genes_info_full.to_csv(Path.joinpath(write_dir, f"{idx}_de_genes_info.csv"))
