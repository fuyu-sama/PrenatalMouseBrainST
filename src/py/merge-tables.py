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
from pathlib import Path

import pandas as pd

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

# %%
for idx in idx_full:
    dfs = {}
    dfs["ai"] = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
        index_col=0,
        header=0,
    ).sort_values(by="Ai", ascending=False)
    all_genes = list(dfs["ai"].index)
    dfs["moran"] = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
        index_col=0,
        header=0,
    )
    dfs["moran"] = dfs["moran"][~dfs["moran"].index.duplicated(keep="last")]
    dfs["moran"] = dfs["moran"].reindex(index=all_genes).sort_values(
        by="I_value", ascending=False)

    dfs["spark"] = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/spark/{idx}-spark.csv"),
        index_col=0,
        header=0,
    )
    dfs["spark"] = dfs["spark"][~dfs["spark"].index.duplicated(keep="last")]
    dfs["spark"] = dfs["spark"].reindex(index=all_genes).sort_values(
        by="adjusted_pvalue", ascending=False)

    dfs["spatialde"] = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/spatialde/{idx}-spatialde.csv"),
        index_col=3,
        header=0,
    )
    dfs["spatialde"] = dfs["spatialde"][~dfs["spatialde"].index.duplicated(
        keep="last")]
    dfs["spatialde"] = dfs["spatialde"].reindex(index=all_genes).sort_values(
        by="qval", ascending=True)
    dfs["somde"] = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/somde/{idx}-somde.csv"),
        index_col=3,
        header=0,
    )
    dfs["somde"] = dfs["somde"][~dfs["somde"].index.duplicated(
        keep="last")]
    dfs["somde"] = dfs["somde"].reindex(index=all_genes).sort_values(
        by="qval", ascending=True)

    merged_df = pd.DataFrame(index=all_genes)
    merged_df["AI"] = dfs["ai"]["Ai"]
    merged_df["AI_rank"] = [dfs["ai"].index.get_loc(i) + 1 for i in all_genes]
    merged_df["global_I"] = dfs["moran"]["I_value"]
    merged_df["I_rank"] = [
        dfs["moran"].index.get_loc(i) + 1 for i in all_genes
    ]
    merged_df["spark_p"] = dfs["spark"]["adjusted_pvalue"]
    merged_df["spark_rank"] = [
        dfs["spark"].index.get_loc(i) + 1 for i in all_genes
    ]
    merged_df["spatialde_q"] = dfs["spatialde"]["qval"]
    merged_df["spatialde_rank"] = [
        dfs["spatialde"].index.get_loc(i) + 1 for i in all_genes
    ]
    merged_df["somde_q"] = dfs["somde"]["qval"]
    merged_df["somde_rank"] = [
        dfs["somde"].index.get_loc(i) + 1 for i in all_genes
    ]
    merged_df = merged_df.sort_values(by="AI_rank")
    merged_df.to_csv(
        Path.joinpath(WORKDIR, f"results/rank-merge/{idx}-merged.csv"),
        na_rep="N/A",
    )
