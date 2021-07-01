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
import gzip
import os
from pathlib import Path

WORKDIR = Path.joinpath(Path.home(), "workspace/mouse-brain-full/")

# %% work
idx = "E165A"

r1_path = Path.joinpath(
    WORKDIR,
    f"Data/{idx}/{idx}_S1_L001_R1_001.fastq.gz.raw.gz",
)
r2_path = Path.joinpath(WORKDIR, f"Data/{idx}/{idx}_S1_L001_R2_001.fastq.gz")

r1_fin = gzip.open(r1_path, "rt")
r2_fin = gzip.open(r2_path, "rt")

os.mkdir(Path.joinpath(WORKDIR, f"Data/{idx}/separate"))
bc_dict = {}
name_list = ["", ""]
seq_list = ["", ""]
id_list = ["", ""]
flag = 1
while 1:
    r1_line = r1_fin.readline()
    r2_line = r2_fin.readline()
    if (not r1_line) or (not r2_line):
        break
    if flag == 1:
        name_list[0] = r1_line
        name_list[1] = r2_line
        flag += 1
    elif flag == 2:
        seq_list[0] = r1_line
        seq_list[1] = r2_line
        flag += 1
    elif flag == 3:
        id_list[0] = r1_line
        id_list[1] = r2_line
        flag += 1
    elif flag == 4:
        bc = seq_list[0][0:16]
        umi = seq_list[0][16:29]
        write_dir = Path.joinpath(WORKDIR, f"Data/{idx}/separate/{bc}")
        if bc not in bc_dict:
            bc_dict[bc] = {}
            os.mkdir(write_dir)
        if umi not in bc_dict[bc]:
            bc_dict[bc][umi] = 1
        r1_fout = gzip.open(
            Path.joinpath(write_dir, f"{umi}_R1.fastq.gz"),
            "at",
        )
        r2_fout = gzip.open(
            Path.joinpath(write_dir, f"{umi}_R2.fastq.gz"),
            "at",
        )
        r1_fout.write(name_list[0] + seq_list[0] + id_list[0] + r1_line)
        r2_fout.write(name_list[1] + seq_list[1] + id_list[1] + r2_line)
        r1_fout.close()
        r2_fout.close()
        # reset variables
        name_list = ["", ""]
        seq_list = ["", ""]
        id_list = ["", ""]
        flag = 1
