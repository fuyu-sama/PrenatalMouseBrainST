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

import sys
from pathlib import Path

fin_path = sys.argv[1]

gene_names = {}
gtf_path = Path.joinpath(
    Path.home(),
    "refgenome/spaceranger/refdata-gex-mm10-2020-A/genes/genes.gtf")
with open(gtf_path) as f:
    for line in f:
        if line[0] == "#":
            continue
        line = line.replace("\n", "")
        line = line.split("\t")
        if line[2] != "gene":
            continue
        line = line[-1].replace("; ", ";")
        line = line.split(";")
        gene_id = ""
        gene_name = ""
        for notes in line:
            note = notes.split(" ")
            if note[0] == "gene_id":
                gene_id = note[1]
            elif note[0] == "gene_name":
                gene_name = note[1]
        if gene_name != "" and gene_id != "":
            gene_name = gene_name.replace('"', "")
            gene_id = gene_id.replace('"', "")
            gene_names[gene_name] = gene_id

with open(fin_path) as fin:
    with open(f"{fin_path}.out", "w") as fout:
        for line in fin:
            if line[0] == ",":
                fout.write("gene_id\n")
                continue
            line = line.replace("\n", "")
            line = line.split(",")
            fout.write(f"{gene_names[line[0]]}\n")
