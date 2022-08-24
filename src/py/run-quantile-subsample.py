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
union_500 = []
union_1000 = []
union_500_99 = []
union_1000_99 = []
for idx in idx_full:
    count_path = Path.joinpath(
        WORKDIR,
        f"Data/scale_df/logcpm/{idx}-logcpm.csv",
    )
    count_df = pd.read_csv(count_path, index_col=0, header=0).T
    mean_series = count_df.mean()
    q = mean_series.quantile(0.99)
    remain_genes = mean_series[mean_series < q].index
    with open(
            Path.joinpath(WORKDIR, f"Data/gene-lists/{idx}-logcpm-0.99.csv"),
            "w",
    ) as f:
        for i in remain_genes:
            f.write("".join([i, "\n"]))

    ai_df = pd.read_csv(
        Path.joinpath(WORKDIR, f"results/Ai/{idx}-Ai.csv"),
        header=0,
        index_col=0,
    ).sort_values(by="Ai", ascending=False)

    for step in [100, 500]:
        for n in range(0, 2500 + step, step):
            selected_genes = ai_df.index[n:n + step]
            with open(
                    Path.joinpath(
                        WORKDIR,
                        f"Data/gene-lists/{idx}-Ai-{n}_{n + step}.csv",
                    ),
                    "w",
            ) as f:
                for i in selected_genes:
                    f.write("".join([i, "\n"]))
    [union_500.append(i) for i in ai_df.index[:500]]
    [union_1000.append(i) for i in ai_df.index[:1000]]

    ai_df = ai_df.reindex(index=[i for i in ai_df.index if i in remain_genes])
    for step in [100, 500]:
        for n in range(0, 2500 + step, step):
            selected_genes = ai_df.index[n:n + step]
            with open(
                    Path.joinpath(
                        WORKDIR,
                        f"Data/gene-lists/{idx}-logcpm-0.99-Ai-{n}_{n + step}.csv",
                    ),
                    "w",
            ) as f:
                for i in selected_genes:
                    f.write("".join([i, "\n"]))
    [union_500_99.append(i) for i in ai_df.index[:500]]
    [union_1000_99.append(i) for i in ai_df.index[:1000]]

with open(
        Path.joinpath(
            WORKDIR,
            f"Data/gene-lists/logcpm-0.99-Ai-500union.csv",
        ),
        "w",
) as f:
    for i in union_500_99:
        f.write("".join([i, "\n"]))
with open(
        Path.joinpath(
            WORKDIR,
            f"Data/gene-lists/logcpm-0.99-Ai-1000union.csv",
        ),
        "w",
) as f:
    for i in union_1000_99:
        f.write("".join([i, "\n"]))
with open(
        Path.joinpath(
            WORKDIR,
            f"Data/gene-lists/Ai-500union.csv",
        ),
        "w",
) as f:
    for i in union_500:
        f.write("".join([i, "\n"]))
with open(
        Path.joinpath(
            WORKDIR,
            f"Data/gene-lists/Ai-1000union.csv",
        ),
        "w",
) as f:
    for i in union_1000:
        f.write("".join([i, "\n"]))
