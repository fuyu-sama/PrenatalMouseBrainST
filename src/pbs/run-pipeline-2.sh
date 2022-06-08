#PBS -N pipeline-2
#PBS -l nodes=comput2:ppn=40
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
scale_methods=(
    combat
)

cd $HOME/workspace/mouse-brain-full

# %% scale data
for directory in ${scale_methods[@]}; do
    if [ ! -d Data/scale_df/${directory} ]; then
        mkdir Data/scale_df/${directory}
    fi
done

if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Scaling data..."
    ${PYTHON_PATH} src/py/run-cpm.py &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-combat.py &>> log/pipeline-2.log
fi

# %% global moran
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Calculating global moran..."
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-global_moran.py \
            ${idx} logcpm 6 &>> log/pipeline-2.log
    done
fi

# %% hotspot
knn=6
scale_method="logcpm"
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Running hotspot..."
    if [ ! -d Data/scale_df/${scale_method}-hotspot-${knn} ]; then
        mkdir Data/scale_df/${scale_method}-hotspot-${knn}
    fi
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-hotspot.py \
            ${idx} ${scale_method} ${knn} &>> log/pipeline-2.log
    done
fi

# %% gene cluster
scale_method="logcpm-hotspot-6"
gene_list=results/gene-cluster/logcpm-hotspot-6-500-union/E135B-11/tables/E135B-genes-11.csv
gene_list_suffix=500-union-E135B-11-11
writedir=results/gene-cluster/${scale_method}-${gene_list_suffix}
if [ ! -d ${writedir} ]; then
    mkdir ${writedir}
fi
if false; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering genes..."
    for idx in ${idx_full[@]}; do
        for n_gene_clusters in {1..10}; do
            if [ ! -d ${writedir}/${idx}-${n_gene_clusters} ]; then
                mkdir ${writedir}/${idx}-${n_gene_clusters}
                mkdir ${writedir}/${idx}-${n_gene_clusters}/tables
            fi
            (
                ${PYTHON_PATH} src/py/run-gene-cluster.py \
                    ${scale_method} ${idx} ${n_gene_clusters} \
                    ${gene_list} ${gene_list_suffix} \
                &>> log/pipeline-2.log
            )&
        done
        wait
    done
fi

# %% cluster
echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with SC3..."
for scale_method in ${scale_methods[@]}; do
    if [ ! -d results/cluster/${scale_method}-sc3 ]; then
        mkdir results/cluster/${scale_method}-sc3
        mkdir results/cluster/${scale_method}-sc3/pattern
        mkdir results/cluster/${scale_method}-sc3/separate
        mkdir results/cluster/${scale_method}-sc3/together
    fi
    for idx in ${idx_full[@]}; do
        (
            Rscript src/R/run-sc3.R ${idx} ${scale_method} \
                &>> log/pipeline-2.log;
            ${PYTHON_PATH} src/py/draw-sc3.py ${idx} ${scale_method} \
                &>> log/pipeline-2.log;
        )&
    done
    wait
done

# %%
finish_time=`date +%Y.%m.%d\ %H:%M:%S`
echo "[${finish_time}] Pipeline-2 finished."
echo "[${finish_time}] Manually check clustering result"
echo "[${finish_time}] and run merge-clusters.py"
