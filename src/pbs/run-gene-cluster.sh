#PBS -N gene-cluster
#PBS -l nodes=comput2:ppn=10
#PBS -l walltime=240:00:00

# %% environment config
PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python
idx_full=(
    E135A E135B 
    E155A E155B 
    E165A E165B 
    E175A1 E175A2 E175B
    P0A1 P0A2
)

cd $HOME/workspace/mouse-brain-full

# %%
scale_method="logcpm-hotspot-6"
gene_list=Data/基因列表1-Thalamus-P0A1-12.csv
gene_list_suffix=基因列表1-Thalamus-P0A1-12
writedir=results/gene-cluster/${scale_method}-${gene_list_suffix}
if [ ! -d ${writedir} ]; then
    mkdir ${writedir}
fi
for idx in ${idx_full[@]}; do
    for n_gene_clusters in {1..20}; do
        if [ ! -d ${writedir}/${idx}-${n_gene_clusters} ]; then
            mkdir ${writedir}/${idx}-${n_gene_clusters}
            mkdir ${writedir}/${idx}-${n_gene_clusters}/tables
        fi
        (
            ${PYTHON_PATH} src/py/run-gene-cluster.py \
                ${scale_method} ${idx} ${n_gene_clusters} \
                ${gene_list} ${gene_list_suffix} \
        )&
    done
    wait
done
