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
import os
import sys
from pathlib import Path

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

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
    "E175B": "V10M17-085-E175B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
timepoints = {
    "E135": ["E135A", "E135B"],
    "E155": ["E155A", "E155B"],
    "E165": ["E165A", "E165B"],
    "E175": ["E175A1", "E175A2", "E175B"],
    "P0": ["P0A1", "P0A2"],
}

regions = ["cortex", "thalamus", "hypothalamus"]

try:
    scale_method = sys.argv[1]
    cluster_method = sys.argv[2]
except IndexError:
    scale_method = "combat-zjq"
    cluster_method = "sc3"

if False:
    for idx in idx_full:
        for region in regions:
            os.mkdir(f"/mnt/Data/mouse-brain-full/spagene/{idx}-{region}")
            read_df = pd.read_csv(
                Path.joinpath(
                    WORKDIR,
                    f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-count.csv",
                ),
                index_col=0,
                header=0,
            )
            for gene in read_df.index:
                os.system(
                    f"cp /mnt/Data/mouse-brain-full/draw_genes/all/{gene}.jpg /mnt/Data/mouse-brain-full/spagene/{idx}-{region}"
                )

if True:
    for tp in timepoints:
        for region in regions:
            for ncs in range(2, 7):
                sets = []
                for idx in timepoints[tp]:
                    read_df = pd.read_csv(
                        Path.joinpath(
                            WORKDIR,
                            f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/tables/{idx}-{region}-{ncs}-de.csv",
                        ),
                        index_col=0,
                        header=0,
                    )
                    read_df = read_df[(read_df["avg_log2FC"] > 0)
                                      & (read_df["p_val_adj"] <= 0.01)]
                    sets.append(set(read_df.index))
                fig, ax = plt.subplots(figsize=(10, 10))
                if len(timepoints[tp]) == 2:
                    venn2(sets, timepoints[tp], ax=ax)
                if len(timepoints[tp]) == 3:
                    venn3(sets, timepoints[tp], ax=ax)
                ax.set_title(f"{tp} {region} {ncs}")
                fig.savefig(
                    Path.joinpath(
                        WORKDIR,
                        f"results/cluster/{scale_method}-{cluster_method}/sub-cluster/{tp}-{region}-{ncs}.jpg",
                    ))
