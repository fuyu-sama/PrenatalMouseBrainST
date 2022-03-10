#PBS -N pipeline-2
#PBS -l nodes=1:ppn=40
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python
idx_full=(E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2)

cd $HOME/workspace/mouse-brain-full

# scale data
if false; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Scaling data..."
    for directory in raw logcpm cpm combat combat-zjq cpm-gmm-2 cpm-gmm-3 combat-gmm-2 combat-gmm-3 logcpm-gmm-2 logcpm-gmm-3; do
        if [ ! -d Data/scale_df/${directory} ]
        then
            mkdir Data/scale_df/${directory}
        fi
    done
    ${PYTHON_PATH} src/py/run-cpm.py &>> log/pipeline-2.log
    ${PYTHON_PATH} src/py/run-combat.py &>> log/pipeline-2.log
    for idx in ${idx_full[@]}; do
        (${PYTHON_PATH} src/py/run-global_moran.py ${idx} cpm 8 &>> log/pipeline-2.log)&
        (${PYTHON_PATH} src/py/run-global_moran.py ${idx} logcpm 8 &>> log/pipeline-2.log)&
        (${PYTHON_PATH} src/py/run-global_moran.py ${idx} combat 8 &>> log/pipeline-2.log)&
        wait
    done
fi

# subsample
# zjq
# ${PYTHON_PATH} src/py/run-subsample.py &>> log/pipeline-2.log
# I gmm
${PYTHON_PATH} src/test-scripts/run-gmm.py cpm &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-gmm.py logcpm &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-gmm.py combat &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py cpm 2 &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py cpm 3 &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py logcpm 2 &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py logcpm 3 &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py combat 2 &>> log/pipeline-2.log
${PYTHON_PATH} src/test-scripts/run-subsample-gmm.py combat 3 &>> log/pipeline-2.log

# cluster
for scale_method in combat-gmm-2 combat-gmm-3 cpm-gmm-2 cpm-gmm-3 logcpm-gmm-2 logcpm-gmm-3; do
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

echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-2 finished."
echo "[`date +%Y.%m.%d\ %H:%M:%S`] Manually check clustering result and run merge-clusters.py"
