#PBS -N pipeline-2
#PBS -l nodes=1:ppn=40
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

# scale data
if false; then
    echo "[`date +%Y.%m.%d\ %H:%M:%S`] Scaling data with combat and Seurat..."
    for directory in raw logcpm combat seurat_integrate z_score; do
        if [ ! -d Data/scale_df/${directory} ]
        then
            mkdir Data/scale_df/${directory}
        fi
    done
    ${PYTHON_PATH} src/py/run-cpm.py &>> log/pipeline-2.log
    (${PYTHON_PATH} src/py/run-combat.py &>> log/pipeline-2.log)&
    (Rscript src/R/run-seurat_integrate.R &>> log/pipeline-2.log)&
    wait
fi

# cluster
for scale_method in combat-qq-logcpm combat-qq-raw; do
    if false; then
        echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with KMeans on all features..."
        if [ ! -d results/cluster/${scale_method}-kmeans ]
        then
            mkdir results/cluster/${scale_method}-kmeans
            mkdir results/cluster/${scale_method}-kmeans/pattern
            mkdir results/cluster/${scale_method}-kmeans/separate
            mkdir results/cluster/${scale_method}-kmeans/together
        fi
        for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
            (${PYTHON_PATH} src/py/run-kmeans.py \
                ${idx} ${scale_method} 1 0 0 &>> log/pipeline-2.log)&
        done
    fi

    if false; then
        echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with KMeans on PCA features..."
        if [ ! -d results/cluster/${scale_method}-pca-kmeans ]; then
            mkdir results/cluster/${scale_method}-pca-kmeans
            mkdir results/cluster/${scale_method}-pca-kmeans/pattern
            mkdir results/cluster/${scale_method}-pca-kmeans/separate
            mkdir results/cluster/${scale_method}-pca-kmeans/together
        fi
        for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
            (${PYTHON_PATH} src/py/run-kmeans.py \
                ${idx} ${scale_method} 0 1 0 &>> log/pipeline-2.log)&
        done
    fi

    if false; then
        echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with KMeans on ICA features..."
        if [ ! -d results/cluster/${scale_method}-ica-kmeans ]; then
            mkdir results/cluster/${scale_method}-ica-kmeans
            mkdir results/cluster/${scale_method}-ica-kmeans/pattern
            mkdir results/cluster/${scale_method}-ica-kmeans/separate
            mkdir results/cluster/${scale_method}-ica-kmeans/together
        fi
        for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
            (
                ${PYTHON_PATH} src/py/run-kmeans.py \
                    ${idx} ${scale_method} 0 0 1 &>> log/pipeline-2.log;
            )&
        done
    fi

    if true; then
        echo "[`date +%Y.%m.%d\ %H:%M:%S`] Clustering with SC3..."
        if [ ! -d results/cluster/${scale_method}-sc3 ]; then
            mkdir results/cluster/${scale_method}-sc3
            mkdir results/cluster/${scale_method}-sc3/pattern
            mkdir results/cluster/${scale_method}-sc3/separate
            mkdir results/cluster/${scale_method}-sc3/together
        fi
        for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2 P0B; do
            (
                Rscript src/R/run-sc3.R ${idx} ${scale_method} &>> log/pipeline-2.log;
                ${PYTHON_PATH} src/py/draw-sc3.py ${idx} ${scale_method} &>> log/pipeline-2.log;
            )&
        done
    fi
done
wait

echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-2 finished."
echo "[`date +%Y.%m.%d\ %H:%M:%S`] Manually check clustering result and run merge-clusters.py"
