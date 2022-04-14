#PBS -N pipeline-2
#PBS -l nodes=comput2:ppn=60
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
    combat
    combat-logcpm-gmm-2
    combat-logcpm-gmm-3 
    combat-logcpm-1000
    logcpm
    cpm
    raw
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

# %% subsample
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Calculating global moran..."
    for scale_method in raw logcpm cpm combat; do
        for idx in ${idx_full[@]}; do
            (${PYTHON_PATH} src/py/run-global_moran.py \
                ${idx} ${scale_method} 8 &>> log/pipeline-2.log)&
        done
        wait
    done
fi

if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling data with I-value..."
    ${PYTHON_PATH} src/py/run-subsample-I.py combat 1000 &>> log/pipeline-2.log
fi

if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling data with I-value & GMM..."
    ${PYTHON_PATH} src/py/run-gmm.py logcpm &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-subsample-gmm.py combat 2 &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-subsample-gmm.py combat 3 &>> log/pipeline-2.log
fi

# %% hotspot
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Running hotspot..."
    for scale_method in raw logcpm cpm combat; do
        if [ ! -d Data/scale_df/${scale_method}-hotspot-8 ]; then
            mkdir Data/scale_df/${scale_method}-hotspot-8
        fi
        for idx in ${idx_full[@]}; do
            (${PYTHON_PATH} src/py/run-hotspot.py \
                ${idx} ${scale_method} 8 &>> log/pipeline-2.log)&
        done
        wait
    done
fi

# %% gene cluster
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering genes..."
    for idx in ${idx_full[@]}; do
        if [ ! -d results/gene-cluster/${idx} ]; then
            mkdir results/gene-cluster/${idx}
        fi
        (${PYTHON_PATH} src/py/run-gene-cluster.py \
            logcpm ${idx} &>> log/pipeline-2.log)
    done
    wait
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
