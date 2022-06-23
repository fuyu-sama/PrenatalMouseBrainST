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
import json
import sys
from pathlib import Path

import pandas as pd
import session_info
from matplotlib import pyplot as plt
from PIL import Image

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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

try:
    idx = sys.argv[1]
    region = sys.argv[2]
    scale_method = sys.argv[3]
except IndexError:
    idx = "E165A"
    region = "cortex"
    scale_method = "combat-Ai-500union"

# %% read data
color_path = Path.joinpath(
    WORKDIR, f"results/cluster/{scale_method}-sc3/color_swift.json")
with open(color_path) as f:
    color_swift = json.load(f)

regions_color = {
    "cortex": ["#CC163A", "#00FFFF"],
    "hypothalamus": ["#EE82EE"],
}
if idx not in color_swift:
    sys.exit()

he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
he_image = Image.open(he_path)
coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
coor_df = pd.read_csv(coor_path, index_col=0, header=0)

cluster_path = Path.joinpath(
    WORKDIR,
    f"results/cluster/{scale_method}-sc3/pattern/{idx}-sc3.csv",
)
cluster_df = pd.read_csv(
    cluster_path,
    index_col=0,
    header=0,
)

in_region = []
for i in cluster_df.index:
    cluster = str(cluster_df.loc[i, "sc3_12_clusters"])
    if color_swift[idx][cluster] in regions_color[region]:
        in_region.append(i)

results_path = Path.joinpath(
    WORKDIR, f"results/RCTD/{region}/{idx}/results/results_df.csv")
results_df = pd.read_csv(results_path, index_col=0, header=0)
weights_path = Path.joinpath(
    WORKDIR, f"results/RCTD/{region}/{idx}/results/weights.csv")
weights_df = pd.read_csv(weights_path, index_col=0, header=0)

coor_df = coor_df.reindex(index=in_region)
results_df = results_df.reindex(index=in_region)
weights_df = weights_df.reindex(index=in_region)

# %% draw weight
for cell_type in weights_df.columns:
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.axis("off")
    ax.imshow(he_image)
    sc = ax.scatter(
        coor_df["X"],
        coor_df["Y"],
        c=weights_df[cell_type],
        cmap="autumn_r",
        s=16,
        alpha=0.7,
    )
    ax.set_title(f"{idx} {cell_type}")
    fig.colorbar(sc, ax=ax)
    fig.savefig(
        Path.joinpath(
            WORKDIR,
            f"results/RCTD/{region}/{idx}/weights/{cell_type}.jpg",
        ))
    plt.close(fig)

# %% draw first type
draw_results = results_df[results_df["spot_class"] != "reject"]
draw_coor = coor_df.reindex(index=draw_results.index)
fig, ax = plt.subplots(figsize=(10, 10))
ax.axis("off")
ax.imshow(he_image)
all_types = list(set(draw_results["first_type"]))
for i in range(len(all_types)):
    draw_results_sub = draw_results[draw_results["first_type"] == all_types[i]]
    draw_coor_sub = draw_coor.reindex(index=draw_results_sub.index)
    ax.scatter(
        draw_coor_sub["X"],
        draw_coor_sub["Y"],
        color=plt.cm.tab20.colors[i],
        label=all_types[i],
        s=16,
    )
ax.set_title(f"{idx} first type")
ax.legend(
    fontsize=16,
    markerscale=2,
    loc="upper left",
    bbox_to_anchor=(1.05, 1),
)
fig.savefig(
    Path.joinpath(WORKDIR, f"results/RCTD/{region}/{idx}/first_type.jpg"),
    bbox_inches="tight",
)
plt.close(fig)

# %% draw second type
draw_results = results_df[results_df["spot_class"] != "reject"]
draw_coor = coor_df.reindex(index=draw_results.index)
fig, ax = plt.subplots(figsize=(10, 10))
ax.axis("off")
ax.imshow(he_image)
all_types = list(set(draw_results["second_type"]))
for i in range(len(all_types)):
    draw_results_sub = draw_results[draw_results["second_type"] ==
                                    all_types[i]]
    draw_coor_sub = draw_coor.reindex(index=draw_results_sub.index)
    ax.scatter(
        draw_coor_sub["X"],
        draw_coor_sub["Y"],
        color=plt.cm.tab20.colors[i],
        label=all_types[i],
        s=16,
    )
ax.set_title(f"{idx} second type")
ax.legend(
    fontsize=16,
    markerscale=2,
    loc="upper left",
    bbox_to_anchor=(1.05, 1),
)
fig.savefig(
    Path.joinpath(WORKDIR, f"results/RCTD/{region}/{idx}/second_type.jpg"),
    bbox_inches="tight",
)
plt.close(fig)
