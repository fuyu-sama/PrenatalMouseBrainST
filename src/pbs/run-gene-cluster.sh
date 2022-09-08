#PBS -N gene-cluster
#PBS -l nodes=comput2:ppn=10
#PBS -l walltime=240:00:00

# %% environment config
PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python
idx_full=(
    E135A E135B
    E155A E155B
    E165A E165B
    E175A1 E175A2
    P0A1 P0A2
)

cd $HOME/workspace/mouse-brain-full

# %%
scale_method="logcpm-hotspot-6"
gene_list=Data/gene-lists/Ai-500union.csv
gene_list_suffix=Ai-500union
if false; then
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
            (${PYTHON_PATH} src/py/run-gene-cluster.py \
                    ${scale_method} ${idx} ${n_gene_clusters} \
                    ${gene_list} ${gene_list_suffix} \
            )&
        done
        wait
    done
    (${PYTHON_PATH} src/py/run-svg_region.py ${scale_method} ${gene_list_suffix})&
    wait
fi

# %%
scale_method="logcpm-hotspot-6"
step_len=500
end_point=2000
if true; then
    for n in `seq 0 ${step_len} ${end_point}`; do
        gene_list=logcpm-0.99-Ai-${n}_`expr ${n} + ${step_len}`
        writedir=results/gene-cluster/${scale_method}-${gene_list}
        if [ ! -d ${writedir} ]; then
            mkdir ${writedir}
        fi
        for idx in ${idx_full[@]}; do
            for n_gene_clusters in {8..20}; do
                if [ ! -d ${writedir}/${idx}-${n_gene_clusters} ]; then
                    mkdir ${writedir}/${idx}-${n_gene_clusters}
                    mkdir ${writedir}/${idx}-${n_gene_clusters}/tables
                fi
                (${PYTHON_PATH} src/py/run-gene-cluster.py \
                    ${scale_method} ${idx} ${n_gene_clusters} \
                    Data/gene-lists/${idx}-${gene_list}.csv ${gene_list})&
            done
            wait
        done
        (${PYTHON_PATH} src/py/run-svg_region.py \
            ${scale_method} ${gene_list} 0.4)&
        (${PYTHON_PATH} src/py/run-svg_region.py \
            ${scale_method} ${gene_list} 0.35)&
        (${PYTHON_PATH} src/py/run-svg_region.py \
            ${scale_method} ${gene_list} 0.3)&
        (${PYTHON_PATH} src/py/run-svg_region.py \
            ${scale_method} ${gene_list} 0.25)&
        (${PYTHON_PATH} src/py/run-svg_region.py \
            ${scale_method} ${gene_list} 0.2)&
        wait
    done
fi
