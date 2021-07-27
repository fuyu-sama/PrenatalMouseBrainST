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

# %% read counts
scale_method = sys.argv[1]
count_path = Path.joinpath(
    WORKDIR,
    f"Data/scale_df/{scale_method}/full-{scale_method}-inter.csv",
)
count_full_df = pd.read_csv(
    count_path,
    index_col=0,
    header=0,
).T

# %% read cluster result
cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/SCT-SC3/pattern/full-SC3.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)
count_full_df = count_full_df.reindex(cluster_df.index)

regions = dict(
    cortex=("E135A_1", "E135B_3", "E135B_9", "E155A_4", "E155B_5", "E155B_6",
            "E165A_1", "E165B_4", "E175A1_8", "E175A2_5", "E175A2_6",
            "E175B_6", "P0A1_5", "P0A2_1", "P0B_3"),
    thalamus=("E135A_2", "E135B_6", "E135B_7", "E155A_5", "E155A_7", "E155A_8",
              "E155B_7", "E165A_3", "E165A_5", "E165B_3", "E175A1_2",
              "E175A1_7", "E175A2_10", "E175A2_11", "E175A2_12", "E175B_2",
              "E175B_4", "E175B_5", "P0A1_6", "P0A1_7", "P0A1_12", "P0A2_6",
              "P0B_1", "P0B_7", "P0B_8"),
    hypothalamus=("E135A_3", "E135A_6", "E135B_5", "E135B_12", "E155A_6",
                  "E155B_1", "E165A_6", "E165B_1", "E175B_3", "E175A1_1",
                  "E175A2_1", "P0A1_1", "P0A2_4", "P0B_12"),
    olfactory=("E135A_5", "E135A_7", "E135B_4", "E155A_11", "E155B_2",
               "E155B_11", "E165A_10", "E165A_11", "E165B_7", "E165B_10",
               "E175A1_5", "E175A2_7", "E175A2_9", "E175B_1", "E175B_9",
               "P0A1_10", "P0A2_7", "P0B_5", "P0B_6"),
    hippocampus=("E155A_2", "E165A_2", "E165B_6", "E175A2_14", "E175A2_15",
                 "E175B_11", "P0A1_8", "P0A2_2", "P0B_4"),
)
regions_label = dict(
    cortex=0,
    thalamus=1,
    hypothalamus=2,
    olfactory=3,
    hippocampus=4,
)
in_regions = [j for i in regions.values() for j in i]
others = [
    i for i in list(set(cluster_df[f"sc3_clusters"])) if i not in in_regions
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
    i for i in umap_df.index if cluster_df.loc[i, "sc3_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, ["red", "green", "blue", "yellow", "orange"]):
    draw_df = umap_df.reindex(index=[
        i for i in umap_df.index
        if cluster_df.loc[i, "sc3_clusters"] in regions[region]
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
        f"results/dimension_reduction/{scale_method}/umap-{scale_method}-1.jpg")
)
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
        f"results/dimension_reduction/{scale_method}/umap-{scale_method}-2.jpg")
)
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
    i for i in densmap_df.index if cluster_df.loc[i, "sc3_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, ["red", "green", "blue", "yellow", "orange"]):
    draw_df = densmap_df.reindex(index=[
        i for i in densmap_df.index
        if cluster_df.loc[i, "sc3_clusters"] in regions[region]
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
        f"results/dimension_reduction/{scale_method}/densmap-{scale_method}-1.jpg"
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
        f"results/dimension_reduction/{scale_method}/densmap-{scale_method}-2.jpg"
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
    i for i in tsne_df.index if cluster_df.loc[i, "sc3_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, ["red", "green", "blue", "yellow", "orange"]):
    draw_df = tsne_df.reindex(index=[
        i for i in tsne_df.index
        if cluster_df.loc[i, "sc3_clusters"] in regions[region]
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
        f"results/dimension_reduction/{scale_method}/tsne-{scale_method}-1.jpg")
)
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
        f"results/dimension_reduction/{scale_method}/tsne-{scale_method}-2.jpg")
)
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
    i for i in densne_df.index if cluster_df.loc[i, "sc3_clusters"] in others
])
ax.scatter(
    draw_df["X"],
    draw_df["Y"],
    label="Others",
    c="grey",
    s=2,
)
for region, c in zip(regions, ["red", "green", "blue", "yellow", "orange"]):
    draw_df = densne_df.reindex(index=[
        i for i in densne_df.index
        if cluster_df.loc[i, "sc3_clusters"] in regions[region]
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
        f"results/dimension_reduction/{scale_method}/densne-{scale_method}-1.jpg"
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
        f"results/dimension_reduction/{scale_method}/densne-{scale_method}-2.jpg"
    ))
plt.close(fig)
