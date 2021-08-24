#PBS -N pipeline-2
#PBS -l nodes=1:ppn=10
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

for scale_method in combat seurat_integrate
do
    for cluster_method in sc3
    do
        if [ ! -d results/cluster/${scale_method}-${cluster_method}/region ]
        then
            mkdir results/cluster/${scale_method}-${cluster_method}/region
        fi

        if [ ! -d results/dimension_reduction/${scale_method}-${cluster_method} ]
        then
            mkdir results/dimension_reduction/${scale_method}-${cluster_method}
        fi

        if [ ! -d results/DE/${scale_method}-${cluster_method} ]
        then
            mkdir results/DE/${scale_method}-${cluster_method}
            mkdir results/DE/${scale_method}-${cluster_method}/region-specific
            mkdir results/DE/${scale_method}-${cluster_method}/timepoint-specific
        fi

        ${PYTHON_PATH} src/py/merge-cluster.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log
        (${PYTHON_PATH} src/py/draw-region.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (Rscript src/R/run-hierarchical.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (Rscript src/R/run-de-timepoint.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
        ${PYTHON_PATH} src/py/draw-de-timepoint.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (Rscript src/R/run-de-region.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
        ${PYTHON_PATH} src/py/draw-de-region.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)
        (${PYTHON_PATH} src/py/run-dim_redu.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        for region in cortex hippocampus hypothalamus thalamus olfactory
        do
        (
            source homer-4.11.sh;
            ${PYTHON_PATH} src/py/run-id-transfer.py \
                results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv \
                &>> log/pipeline-3.log;

            findMotifs.pl \
                results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv.out \
                mouse results/motifResults/${scale_method}-${cluster_method}/${region} \
                -start -1000 -end 1000 &>> log/pipeline-3.log
        )&
        done
        wait
        (${PYTHON_PATH} src/py/draw-homer.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (${PYTHON_PATH} src/py/draw-go.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done
wait
