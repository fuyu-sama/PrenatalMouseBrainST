#PBS -N pipeline-2
#PBS -l nodes=comput2:ppn=40
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
scale_methods=(
    combat-gmm
    #combat
    #logcpm
    #cpm
    #raw
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
            ${idx} logcpm 8 &>> log/pipeline-2.log
    done
fi

# %% subsample
knn=8
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling data with I-value & GMM..."
    if [ ! -d results/I-gmm/logcpm-${knn} ]; then
        mkdir results/I-gmm/logcpm-${knn}
        mkdir results/I-gmm/logcpm-${knn}/venn
        mkdir results/I-gmm/logcpm-${knn}/pdf
    fi
    ${PYTHON_PATH} src/py/run-gmm.py logcpm ${knn}&>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-subsample-gmm.py combat ${knn} &>> log/pipeline-2.log
fi

# %% hotspot
knn=8
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Running hotspot..."
    if [ ! -d Data/scale_df/logcpm-hotspot-${knn} ]; then
        mkdir Data/scale_df/logcpm-hotspot-${knn}
    fi
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-hotspot.py \
            ${idx} logcpm ${knn} &>> log/pipeline-2.log
    done
fi

# %% gene cluster
knn=8
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering genes..."
    for idx in ${idx_full[@]}; do
        for n_gene_clusters in {6..12}; do
            if [ ! -d results/gene-cluster/${idx}-${n_gene_clusters} ]; then
                mkdir results/gene-cluster/${idx}-${n_gene_clusters}
                mkdir results/gene-cluster/${idx}-${n_gene_clusters}/tables
            fi
            (
                ${PYTHON_PATH} src/py/run-gene-cluster.py \
                    logcpm ${idx} ${knn} ${n_gene_clusters} &>> log/pipeline-2.log
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
