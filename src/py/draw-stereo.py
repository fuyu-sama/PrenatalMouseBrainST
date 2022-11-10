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
from pathlib import Path

import stereo
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans

plt.rcParams.update({"font.size": 64})
WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")

# %%
bin_size = 1
ste = stereo.io.read_gef(
    str(
        Path.joinpath(
            WORKDIR,
            "SAW.4.1.0/SS200000977TL_A1_result/SS200000977TL_A1.gef",
        )),
    bin_size=bin_size,
)

count_df = ste.to_df()
coor_df = pd.DataFrame(ste.position, index=ste.cell_names, columns=["X", "Y"])

# %%
if False:
    model_kmeans = KMeans(n_clusters=3)
    result_kmeans = model_kmeans.fit_predict(
        count_df.T.sum().to_numpy().reshape(-1, 1), )
    result_kmeans = pd.Series(result_kmeans, index=count_df.index)
    max_cluster, max_value = -1, -1
    for i in range(3):
        cluster_length = result_kmeans[result_kmeans == i].shape[0]
        if cluster_length > max_value:
            max_cluster, max_value = i, cluster_length
    count_df = count_df.drop(
        index=result_kmeans[result_kmeans == max_cluster].index)

# %%
genes = [
    "9130024F11Rik", "Mef2c", "Satb2", "Sla", "Sox5", "3110039M20RiK",
    "Bcl11b", "Bhlhe22", "Fezf2", "Foxg1", "Hivep2", "Mpped1", "Neurod2",
    "Neurod6", "Tbr1", "Zbtb18", "Zeb2", "Eomes", "Fbln2", "Gas1", "Gm11266",
    "Insm1", "Mki67", "Neurod1", "Neurog2", "Rlbp1", "Sema3c", "Sstr2", "Tac2",
    "Thrsp", "Top2a", "Zbtb20"
]
genes = ["Neurod2", "Zbtb18"]
write_dir = Path.joinpath(WORKDIR, f"draw_genes/stereo")
if not write_dir.exists():
    write_dir.mkdir()
for gene in genes:
    if gene not in count_df.columns:
        print(gene)
        continue
    fig, ax = plt.subplots(figsize=(43, 40))
    exp_spots = count_df[gene][count_df[gene] > 0].index
    sc = ax.scatter(
        coor_df["X"].reindex(index=exp_spots),
        coor_df["Y"].reindex(index=exp_spots),
        c=count_df[gene].reindex(index=exp_spots),
        s=1,
        cmap="autumn",
        vmin=0,
    )
    ax.set_title(f"{gene}, bin{bin_size}")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(sc)
    fig.savefig(Path.joinpath(write_dir, f"{gene}.jpg"))
