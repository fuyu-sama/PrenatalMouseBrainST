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

from pathlib import Path

import pandas as pd
import session_info

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
    "E175A2": "V10M17-101-E175B2",
    "E175B": "V10M17-085-E175B",
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}

for idx in idx_full:
    print(idx)
    position_path = Path.joinpath(
        WORKDIR, f"spaceranger/{idx}/outs/spatial/tissue_positions_list.csv")
    count_path = Path.joinpath(
        WORKDIR,
        f"spaceranger/{idx}/outs/filtered_feature_bc_matrix/{idx}.csv")
    feature_path = Path.joinpath(
        WORKDIR,
        f"paceranger/{idx}/outs/filtered_feature_bc_matrix/features.tsv.gz")
    position_df = pd.read_csv(position_path, index_col=0, header=None)
    count_df = pd.read_csv(count_path, index_col=0, header=0)
    position_df = position_df.reindex(index=count_df.columns)
    feature_df = pd.read_csv(
        feature_path,
        index_col=0,
        header=None,
        compression="gzip",
        sep="\t",
    )

    genes = []
    for i in count_df.index:
        gene_name = feature_df.loc[i, 1]
        while 1:
            if gene_name in genes:
                gene_name += "*"
            else:
                break
        genes.append(gene_name)
    count_df.index = genes

    coor_df = pd.DataFrame()
    coor_df["X"] = position_df[5]
    coor_df["Y"] = position_df[4]
    count_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"spaceranger/{idx}/outs/filtered_feature_bc_matrix/{idx}.csv",
        ))
    coor_df.to_csv(
        Path.joinpath(
            WORKDIR,
            f"spaceranger/{idx}/outs/spatial/coor-{idx}.csv",
        ))
