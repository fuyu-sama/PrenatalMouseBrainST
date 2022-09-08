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
import math
import sys
from copy import deepcopy
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image
from seaborn import color_palette

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
colors_short = [
    "tab:red", "tab:green", "tab:blue", "tab:orange", "tab:pink", "tab:cyan",
    "k"
]

try:
    scale_method = sys.argv[1]
    suffix = sys.argv[2]
    threshold = float(sys.argv[3])
except IndexError:
    scale_method = "logcpm-hotspot-6"
    suffix = "logcpm-0.99-Ai-0_500"
    threshold = 0.4


# %%
def main(idx):
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

    for ncs in range(8, 21):
        gene_result = pd.read_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-genes.csv",
            ),
            index_col=0,
            header=0,
        )["0"]

        # build mean_df and ai_rank_df
        ai_rank_df = pd.DataFrame()
        mean_df = pd.DataFrame(index=count_df.index)
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
        # END: for i in range(1, ncs + 1):
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

        # build type_df
        columns = ["spot_type"]
        [columns.append(f"type_{i}") for i in range(1, 1 + ncs)]
        type_df = pd.DataFrame(columns=columns)
        for spot in mean_df.index:
            temp_df = pd.DataFrame(index=[spot], columns=columns)
            sort_series = mean_df.loc[spot, :].sort_values(ascending=False)
            n_multi = len(sort_series[sort_series > threshold])
            spot_type = "singlet"
            if sort_series.iloc[0] < 0.2:
                spot_type = "uncertain"
            if n_multi > 1:
                spot_type = f"{n_multi}_multi_types"
            temp_df.loc[spot, "spot_type"] = spot_type
            temp_df.loc[spot, columns[1:]] = sort_series.index
            type_df = pd.concat([type_df, temp_df])
        # END: for spot in mean_df.index:
        type_df.to_csv(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-spots.csv",
            ))

        # draw spots without multi-clusters
        draw_certain = False
        if draw_certain:
            certain_spots = type_df[type_df["spot_type"] != "uncertain"].index
        else:
            certain_spots = type_df.index
        fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
        ax.scatter(
            coor_df["X"].reindex(index=certain_spots),
            coor_df["Y"].reindex(index=certain_spots),
            s=16,
            c=type_df[f"type_1"].reindex(index=certain_spots),
            cmap=ListedColormap(colors),
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
                f"{idx}-{ncs}/{idx}-spots-1.jpg",
            ))
        plt.close(fig)

        # draw spots separate
        nrows = math.ceil(ncs / 5)
        fig, axes = plt.subplots(nrows, 5, figsize=(50, nrows * 10))
        [ax.axis("off") for ax in axes.flatten()]
        for n, ax in zip(range(1, ncs + 1), axes.flatten()):
            in_spots = type_df[type_df[f"type_1"] == n].index
            ax.scatter(
                coor_df["X"].reindex(index=certain_spots).reindex(
                    index=in_spots),
                coor_df["Y"].reindex(index=certain_spots).reindex(
                    index=in_spots),
                s=16,
            )
            ax.imshow(he_image)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{idx} ncs {ncs} cluster {n}")
        # END: for n, ax in zip(range(1, ncs + 1), axes.flatten()):
        fig.savefig(
            Path.joinpath(
                WORKDIR,
                f"results/gene-cluster/{scale_method}-{suffix}/",
                f"{idx}-{ncs}/{idx}-spots-separate.jpg",
            ))
        plt.close(fig)

        # render type_df to dict
        # key: spot, value: set of clusters
        singlet_df = type_df[type_df["spot_type"] == "singlet"]
        spot_type_dict = {}
        for n in range(2, ncs + 1):
            multi_spots = type_df[type_df["spot_type"] == f"{n}_multi_types"]
            for spot in multi_spots.index:
                spot_type_dict[spot] = frozenset([
                    multi_spots.loc[spot, f"type_{i}"]
                    for i in range(1, n + 1)
                ])
        # END: for n in range(2, ncs + 1):

        for n in [2, 3]:
            all_n_combinations = []
            for sets in frozenset(spot_type_dict.values()):
                if len(sets) != n:
                    continue

                # draw multi-cluster map
                fig, ax = plt.subplots(figsize=(10, 10))
                ax.imshow(he_image)
                ax.axis("off")

                all_combinations = []
                n_combinations = 2
                while n_combinations <= n:
                    for c in combinations(sets, n_combinations):
                        c = frozenset(c)
                        all_combinations.append(c)
                        all_n_combinations.append(c)
                    n_combinations += 1

                # single cluster spots
                for i, type_ in enumerate(sets):
                    type_spots = list(
                        singlet_df[singlet_df["type_1"] == type_].index)
                    for spot in spot_type_dict:
                        if type_ in spot_type_dict[spot]:
                            if spot_type_dict[spot] not in all_combinations:
                                type_spots.append(spot)
                    ax.scatter(
                        coor_df["X"].reindex(index=type_spots),
                        coor_df["Y"].reindex(index=type_spots),
                        s=16,
                        c=colors_short[i],
                        label=type_,
                    )
                # END: for i, type_ in enumerate(sets):
                for c in all_combinations:
                    i += 1
                    combination_spots = []
                    for spot in spot_type_dict:
                        if spot_type_dict[spot].issuperset(c):
                            combination_spots.append(spot)
                    ax.scatter(
                        coor_df["X"].reindex(index=combination_spots),
                        coor_df["Y"].reindex(index=combination_spots),
                        s=16,
                        c=colors_short[i],
                        label=" & ".join([str(i) for i in c]),
                    )
                # END: for c in all_combinations:
                fig.legend()
                fig.savefig(
                    Path.joinpath(
                        WORKDIR,
                        f"results/gene-cluster/{scale_method}-{suffix}/{idx}-{ncs}",
                        f"{idx}-spots-{threshold}-{'_'.join(str(i) for i in sets)}.jpg",
                    ))
                plt.close(fig)
            # END: for sets in frozenset(spot_type_dict.values()):

            # draw spots with multi-clusters
            all_n_combinations = frozenset(all_n_combinations)
            label_spot_dict = {}
            for c in all_n_combinations:
                combination_spots = []
                for spot in spot_type_dict:
                    if spot_type_dict[spot].issuperset(c):
                        combination_spots.append(spot)
                if len(combination_spots) > 10:
                    label = " & ".join([str(i) for i in c])
                    label_spot_dict[label] = deepcopy(combination_spots)
            fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
            ax.imshow(he_image)
            ax.axis("off")
            for i, cluster in enumerate(range(1, ncs + 1)):
                region_spots = type_df[type_df["type_1"] == cluster].index
                ax.scatter(
                    coor_df["X"].reindex(index=region_spots),
                    coor_df["Y"].reindex(index=region_spots),
                    s=16,
                    color=colors[i],
                    label=cluster,
                )
            # END: for i, cluster in enumerate(range(1, ncs + 1)):
            auto_colors = color_palette("hls", len(label_spot_dict))
            for i, (label, combination_spots) in enumerate(
                    label_spot_dict.items()):
                ax.scatter(
                    coor_df["X"].reindex(index=combination_spots),
                    coor_df["Y"].reindex(index=combination_spots),
                    s=16,
                    color=auto_colors[i],
                    label=label,
                )
            # END: for i, (label, combination_spots) in enumerate(
            #              label_spot_dict.items()):
            fig.legend(ncol=7)
            fig.savefig(
                Path.joinpath(
                    WORKDIR,
                    f"results/gene-cluster/{scale_method}-{suffix}/",
                    f"{idx}-{ncs}/{idx}-spots-{threshold}-{n}_combinations.jpg",
                ))
            plt.close(fig)
        # END: for n in [2, 3]:
    he_image.close()
    # END: for ncs in range(8, 21):


# END: def main():

# %%
wanted_genes = [
    "Mef2c", "Bcl11a", "Sox5", "Hivep2", "Satb2", "Fezf2", "Neurod2",
    "Neurod6", "Zbtb18", "Tbr1", "Zfp423", "Tcf7l2", "Zic1", "Zic3", "Zic4",
    "Zic5", "Gbx2", "Id4", "Slc18a2", "Clybl", "Lhx9", "Foxp2", "Dlx1", "Otp",
    "Asb4", "Trh", "Ddc", "Dlk1", "Nkx2-2", "Magel2", "Nap1l5", "Peg10"
]
pool = Pool(10)
pool.map(main, idx_full)
pool.close()
pool.join()
