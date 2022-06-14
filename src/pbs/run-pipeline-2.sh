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
    combat-Ai-500union
    combat-Ai-1000union
    combat-Ai-0_500
    combat-Ai-500_1000
    combat-Ai-1000_1500
    combat-Ai-1500_2000
    combat-Ai-2000_2500
    combat-Ai-0_1000
    combat-Ai-0_1500
    combat-Ai-0_2000
    combat-Ai-0_2500
    combat-I-500union
    combat-I-1000union
    combat-I-0_500
    combat-I-500_1000
    combat-I-1000_1500
    combat-I-1500_2000
    combat-I-2000_2500
    combat-I-0_1000
    combat-I-0_1500
    combat-I-0_2000
    combat-I-0_2500
    combat
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

# %% global moran
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Calculating global moran..."
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-global_moran.py \
            ${idx} logcpm 6 &>> log/pipeline-2.log
    done
fi

# %% hotspot
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Running hotspot..."
    if [ ! -d Data/scale_df/logcpm-hotspot-6 ]; then
        mkdir Data/scale_df/logcpm-hotspot-6
    fi
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-hotspot.py \
            ${idx} logcpm 6 &>> log/pipeline-2.log
    done
fi

# %% density
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Running density..."
    for idx in ${idx_full[@]}; do
        ${PYTHON_PATH} src/py/run-density.py ${idx} &>> log/pipeline-2.log
    done
fi

# %% subsample
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling..."
    ${PYTHON_PATH} src/py/run-subsample.py combat &>> log/pipeline-2.log
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
