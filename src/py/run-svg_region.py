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

# %%
import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image

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
    # "E175B": "V10M17-085-E175B",
    # "P0B": "V10M17-100-P0B",
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
    scale_method = sys.argv[1]
    suffix = sys.argv[2]
except IndexError:
    scale_method = "logcpm-hotspot-6"
    suffix = "Ai-500union"

# %%
wanted_genes = [
    "Mef2c", "Bcl11a", "Sox5", "Hivep2", "Satb2", "Fezf2", "Neurod2",
    "Neurod6", "Zbtb18", "Tbr1", "Zfp423", "Tcf7l2", "Zic1", "Zic3", "Zic4",
    "Zic5", "Gbx2", "Id4", "Slc18a2", "Clybl", "Lhx9", "Foxp2", "Dlx1", "Otp",
    "Asb4", "Trh", "Ddc", "Dlk1", "Nkx2-2", "Magel2", "Nap1l5", "Peg10"
]
for idx in idx_full:
    count_df = pd.read_csv(
        Path.joinpath(
            WORKDIR,
            f"Data/scale_df/{scale_method}/{idx}-{scale_method}.csv",
        ),
        index_col=0,
        header=0,
    ).T

    he_path = Path.joinpath(WORKDIR, f"Data/HE/{idx_full[idx]}.tif")
    he_image = Image.open(he_path)
    coor_path = Path.joinpath(WORKDIR, f"Data/coor_df/{idx}-coor.csv")
    coor_df = pd.read_csv(coor_path, index_col=0, header=0)
    assert all(count_df.index == coor_df.index)

    ai_df = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
        index_col=0,
        header=0,
    ).sort_values(by="Ai", ascending=False)

    mean_df = pd.DataFrame(index=count_df.index)
    for ncs in range(2, 21):
        gene_result = pd.read_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-genes.csv",
            ),
            index_col=0,
            header=0,
        )["0"]
        ai_rank_df = pd.DataFrame()
        for i in range(1, ncs + 1):
            genes = gene_result[gene_result == i].index
            mean_df[i] = count_df[genes].T.mean()
            selected_genes = [i for i in genes if i in wanted_genes]
            selected_genes = [i for i in selected_genes if i in ai_df.index]
            if len(selected_genes) > 0:
                write_df = pd.DataFrame()
                write_df["region"] = [i] * len(selected_genes)
                write_df["gene"] = selected_genes
                write_df["AI"] = ai_df.reindex(
                    index=selected_genes,
                    columns=["Ai"],
                ).values
                write_df["rank"] = [
                    ai_df.index.get_loc(i) + 1 for i in selected_genes
                ]
                write_df = write_df.sort_values(by="rank")
                ai_rank_df = pd.concat([ai_rank_df, write_df])
        ai_rank_df.index = [i for i in range(ai_rank_df.shape[0])]
        ai_rank_df.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-rank.csv",
            ))
        mean_df.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-mean.csv",
            ))
        mean_df.T.idxmax().to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-spots.csv",
            ))
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.scatter(
            coor_df["X"],
            coor_df["Y"],
            s=16,
            c=mean_df.T.idxmax(),
            cmap=ListedColormap(colors[:ncs]),
            alpha=0.7,
        )
        ax.imshow(he_image)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{idx} {ncs} clusters")
        ax.axis("off")
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-spots.jpg",
            ))
        plt.close(fig)
    he_image.close()
