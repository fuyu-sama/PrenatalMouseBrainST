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
import sys
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from sklearn.manifold import TSNE
from umap import UMAP

import densne

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")
session_info.show()
plt.rcParams.update({"font.size": 24})

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

scale_method = sys.argv[1]
cluster_method = sys.argv[2]

# %% read counts
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-{cluster_method}/pattern/full-{cluster_method}.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
count_full_df = count_full_df.reindex(cluster_df.index)

regions_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-{cluster_method}/regions.json")
with open(regions_path) as f:
    regions = json.load(f)["regions"]
regions_label = dict(
    cortex=0,
    hippocampus=1,
    thalamus=2,
    hypothalamus=3,
    striatum=4,
    mge=5,
    amygdalar=6,
)
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"{cluster_method}_clusters"]))
    if i not in in_regions
]

# %% umap
model_umap = UMAP(random_state=42)
umap_df = pd.DataFrame(
    model_umap.fit_transform(count_full_df),
    index=count_full_df.index,
    columns=["X", "Y"],
)

# %% color by region
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
draw_df = umap_df.reindex(index=[
    i for i in umap_df.index
    if cluster_df.loc[i, f"{cluster_method}_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, colors):
    draw_df = umap_df.reindex(index=[
        i for i in umap_df.index
        if cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region]
    ])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=region,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10)
ax.set_title(f"UMAP {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/umap-{scale_method}-1.jpg"
    ))
plt.close(fig)

# %% color by sample
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
for tp, c in zip(idx_full.keys(), colors):
    draw_df = umap_df.reindex(index=[i for i in umap_df.index if tp in i])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=tp,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10, ncol=2)
ax.set_title(f"UMAP {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/umap-{scale_method}-2.jpg"
    ))
plt.close(fig)

# %% densmap
model_densmap = UMAP(random_state=42, densmap=True)
densmap_df = pd.DataFrame(
    model_densmap.fit_transform(count_full_df),
    index=count_full_df.index,
    columns=["X", "Y"],
)

# %% color by region
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
draw_df = densmap_df.reindex(index=[
    i for i in densmap_df.index
    if cluster_df.loc[i, f"{cluster_method}_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, colors):
    draw_df = densmap_df.reindex(index=[
        i for i in densmap_df.index
        if cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region]
    ])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=region,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10)
ax.set_title(f"densmap {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/densmap-{scale_method}-1.jpg"
    ))
plt.close(fig)

# %% color by sample
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
for tp, c in zip(idx_full.keys(), colors):
    draw_df = densmap_df.reindex(
        index=[i for i in densmap_df.index if tp in i])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=tp,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10, ncol=2)
ax.set_title(f"densmap {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/densmap-{scale_method}-2.jpg"
    ))
plt.close(fig)

# %% tsne
model_tsne = TSNE(random_state=42)
tsne_df = pd.DataFrame(
    model_tsne.fit_transform(count_full_df),
    index=count_full_df.index,
    columns=["X", "Y"],
)

# %% color by region
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
draw_df = tsne_df.reindex(index=[
    i for i in tsne_df.index
    if cluster_df.loc[i, f"{cluster_method}_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, colors):
    draw_df = tsne_df.reindex(index=[
        i for i in tsne_df.index
        if cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region]
    ])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=region,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10)
ax.set_title(f"tSNE {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/tsne-{scale_method}-1.jpg"
    ))
plt.close(fig)

# %% color by sample
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
for tp, c in zip(idx_full.keys(), colors):
    draw_df = tsne_df.reindex(index=[i for i in tsne_df.index if tp in i])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=tp,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10, ncol=2)
ax.set_title(f"tSNE {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/tsne-{scale_method}-2.jpg"
    ))
plt.close(fig)

# %% densne
result_densne = densne.run_densne(count_full_df)[0]
densne_df = pd.DataFrame(result_densne,
                         index=count_full_df.index,
                         columns=["X", "Y"])

# %% color by region
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
draw_df = densne_df.reindex(index=[
    i for i in densne_df.index
    if cluster_df.loc[i, f"{cluster_method}_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, colors):
    draw_df = densne_df.reindex(index=[
        i for i in densne_df.index
        if cluster_df.loc[i, f"{cluster_method}_clusters"] in regions[region]
    ])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=region,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10)
ax.set_title(f"densne {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/densne-{scale_method}-1.jpg"
    ))
plt.close(fig)

# %% color by sample
fig, ax = plt.subplots(figsize=(20, 20))
ax.set_xticks([])
ax.set_yticks([])
for tp, c in zip(idx_full.keys(), colors):
    draw_df = densne_df.reindex(index=[i for i in densne_df.index if tp in i])
    ax.scatter(
        draw_df["X"],
        draw_df["Y"],
        label=tp,
        c=c,
        s=2,
    )
ax.legend(loc="best", fontsize=16, markerscale=10, ncol=2)
ax.set_title(f"densne {scale_method}")
fig.savefig(
    Path.joinpath(
        WORKDIR,
        f"results/dimension_reduction/{scale_method}-{cluster_method}/densne-{scale_method}-2.jpg"
    ))
plt.close(fig)
