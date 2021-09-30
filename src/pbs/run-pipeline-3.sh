#PBS -N pipeline-2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

# RCTD
for region in hypothalamus cortex; do
    if [ ! -d results/RCTD/${region} ]; then
        mkdir results/RCTD/${region}
    fi
    for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
        (
            if [ ! -d results/RCTD/${region}/${idx} ]; then
                mkdir results/RCTD/${region}/${idx}
                mkdir results/RCTD/${region}/${idx}/results
                mkdir results/RCTD/${region}/${idx}/weights
            fi

            source hdf5-1.12.0.sh
            Rscript src/R/run-rctd-${region}.R \
                ${idx} combat sc3 &>> log/pipeline-3.log;

            exit_status=$?
            if [ $exit_status -eq 0 ]; then
                ${PYTHON_PATH} src/py/draw-rctd.py \
                    ${idx} ${region} & >> log/pipeline-3.log;
            fi
        )&
    done
done

for scale_method in combat; do
    for cluster_method in sc3; do
        if [ ! -d results/cluster/${scale_method}-${cluster_method}/region ]; then
            mkdir results/cluster/${scale_method}-${cluster_method}/region
        fi

        if [ ! -d results/cluster/${scale_method}-${cluster_method}/sub-clusters ]; then
            mkdir results/cluster/${scale_method}-${cluster_method}/sub-clusters
            mkdir results/cluster/${scale_method}-${cluster_method}/sub-clusters/pattern
            mkdir results/cluster/${scale_method}-${cluster_method}/sub-clusters/together
        fi

        if [ ! -d results/dimension_reduction/${scale_method}-${cluster_method} ]; then
            mkdir results/dimension_reduction/${scale_method}-${cluster_method}
        fi

        if [ ! -d results/DE/${scale_method}-${cluster_method} ]; then
            mkdir results/DE/${scale_method}-${cluster_method}
            mkdir results/DE/${scale_method}-${cluster_method}/region-specific
            mkdir results/DE/${scale_method}-${cluster_method}/timepoint-specific
        fi

        # draw regions
        (${PYTHON_PATH} src/py/draw-region.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&

        # dimension-reduction
        (${PYTHON_PATH} src/py/run-dim_redu.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&

        # draw hierarchical clustering tree
        (Rscript src/R/run-hierarchical.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&

        # DE
        (Rscript src/R/run-de-region.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
        ${PYTHON_PATH} src/py/draw-de-region.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
        Rscript src/R/run-de-timepoint.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
        ${PYTHON_PATH} src/py/draw-de-timepoint.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&

        # sub-clustering
        for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
            for region in cortex hypothalamus; do
            (
                Rscript src/R/run-sub-sc3.R \
                    ${idx} ${scale_method} ${region} &>> log/pipeline-3.log;
                ${PYTHON_PATH} src/py/draw-sub-sc3.py \
                    ${idx} ${scale_method} ${region} &>> log/pipeline-3.log;
            )&
            done
        done

        # homer
        for region in cortex hippocampus hypothalamus thalamus amygdalar mge striatum; do
        (
            source homer-4.11.sh;
            ${PYTHON_PATH} src/py/run-id-transfer.py \
                results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv \
                &>> log/pipeline-3.log;

            findMotifs.pl \
                results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv.out \
                mouse results/motifResults/${scale_method}-${cluster_method}/${region} \
                -start -1000 -end 1000 &>> log/pipeline-3.log
        )& pids+=($!)
        done
        wait "${pids[@]}"
        (${PYTHON_PATH} src/py/draw-homer.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (${PYTHON_PATH} src/py/draw-go.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done
wait
