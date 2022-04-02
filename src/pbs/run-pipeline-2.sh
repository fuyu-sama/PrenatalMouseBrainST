#PBS -N pipeline-2
#PBS -l nodes=1:ppn=40
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python
idx_full=(
    E135A E135B 
    E155A E155B 
    E165A E165B 
    E175A1 E175A2 E175B 
    P0A1 P0A2
)
scale_methods=(
    #raw 
    #cpm cpm-gmm-2 cpm-gmm-3 
    #logcpm logcpm-gmm-2 logcpm-gmm-3 
    #combat combat-gmm-2 combat-gmm-3 
    combat-logcpm-1000
)

cd $HOME/workspace/mouse-brain-full

# scale data
for directory in ${scale_methods[@]}; do
    if [ ! -d Data/scale_df/${directory} ]
    then
        mkdir Data/scale_df/${directory}
    fi
done

if false; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Scaling data..."
    ${PYTHON_PATH} src/py/run-cpm.py &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-combat.py &>> log/pipeline-2.log
fi

# subsample
if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Calculating global moran..."
    for idx in ${idx_full[@]}; do
        for scale_method in cpm logcpm combat; do
            ${PYTHON_PATH} src/py/run-global_moran.py \
                ${idx} ${scale_method} 8 &>> log/pipeline-2.log
        done
    done
fi

if true; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling data with I-value..."
    ${PYTHON_PATH} src/py/run-subsample-I.py combat 1000 &>> log/pipeline-2.log
fi

if false; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Subsampling data with I-value & GMM..."
    ${PYTHON_PATH} src/py/run-gmm.py logcpm &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-subsample-gmm.py combat 2 &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-subsample-gmm.py combat 3 &>> log/pipeline-2.log
fi

# cluster
for scale_method in ${scale_methods[@]}; do
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with SC3..."
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
done
wait

finish_time=`date +%Y.%m.%d\ %H:%M:%S`
echo "[${finish_time}] Pipeline-2 finished."
echo "[${finish_time}] Manually check clustering result"
echo "[${finish_time}] and run merge-clusters.py"
