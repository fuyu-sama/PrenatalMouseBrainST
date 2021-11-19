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
from multiprocessing import Pool
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

idx_group = [
    "E135A", "E155A", "E165A", "E175A1", "E175B", "P0A1", "E135B", "E155B",
    "E165B", "E175A2", "P0B", "P0A2"
]

# %% read data
count_dict = {}
coor_dict = {}
he_dict = {}

for idx in idx_group:
    count_path = Path.joinpath(
        WORKDIR,
        f"spaceranger/{idx}/outs/filtered_feature_bc_matrix/{idx}.csv",
    )
    coor_path = Path.joinpath(
        WORKDIR,
        f"spaceranger/{idx}/outs/spatial/coor-{idx}.csv",
    )
    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")

    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    he_image = Image.open(he_path)
    assert all(count_df.index == coor_df.index)

    count_df.index = [f"{idx}_{i}" for i in count_df.index]
    coor_df.index = [f"{idx}_{i}" for i in coor_df.index]

    count_dict[idx] = count_df
    coor_dict[idx] = coor_df
    he_dict[idx] = he_image
genes = count_dict["E135A"].columns


# %% draw
def draw_genes(range_list: list) -> None:
    for i in range(range_list[0], range_list[1]):
        gene = genes[i]
        fig, axes = plt.subplots(2, 6, figsize=(60, 20))
        [ax.set_xticks([]) for ax in axes.flatten()]
        [ax.set_yticks([]) for ax in axes.flatten()]
        [
            axis.set_visible(False) for ax in axes.flatten()
            for axis in ax.spines.values()
        ]
        for idx, ax in zip(idx_group, axes.flatten()):
            draw_counts = count_dict[idx][gene].loc[count_dict[idx][gene] > 0]
            draw_coor = coor_dict[idx].reindex(draw_counts.index)
            ax.imshow(he_dict[idx])
            ax.set_title(f"{idx} {gene}")
            sc = ax.scatter(
                draw_coor["X"],
                draw_coor["Y"],
                c=draw_counts,
                cmap="autumn_r",
                s=16,
                alpha=0.7,
                vmax=count_df[gene].max() + 1,
                vmin=0,
            )
            fig.colorbar(sc, ax=ax)
        fig.savefig(Path.joinpath(WORKDIR, f"draw_genes/all-raw/{gene}.jpg"))
        plt.close(fig)


pool = Pool(5)
pool.map(draw_genes, [[i * 6457, (i + 1) * 6457 + 1] for i in range(5)])
pool.close()
pool.join()
