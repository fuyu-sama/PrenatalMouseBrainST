#PBS -N pipeline-3
#PBS -l nodes=1:ppn=20
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
    raw 
    cpm cpm-gmm-2 cpm-gmm-3 
    logcpm logcpm-gmm-2 logcpm-gmm-3 
    combat combat-gmm-2 combat-gmm-3 
    combat-1000
)
regions=(cortex hippocampus hypothalamus thalamus "amygdalar.olfactory")

cd $HOME/workspace/mouse-brain-full

# RCTD
for region in hypothalamus cortex; do
    if [ ! -d results/RCTD/${region} ]; then
        mkdir results/RCTD/${region}
    fi
    for idx in ${idx_full[@]}; do
        (
            if [ ! -d results/RCTD/${region}/${idx} ]; then
                mkdir results/RCTD/${region}/${idx}
                mkdir results/RCTD/${region}/${idx}/results
                mkdir results/RCTD/${region}/${idx}/weights
            fi

            source hdf5-1.12.0.sh
            Rscript src/R/run-rctd-${region}.R \
                ${idx} combat-gmm-2 sc3 &>> log/pipeline-3.log;

            ${PYTHON_PATH} src/py/draw-rctd.py \
                ${idx} ${region} combat-gmm-2 sc3 &>> log/pipeline-3.log;
            #${PYTHON_PATH} src/py/draw-rctd-IE.py \
                #${idx} ${region} combat-gmm-2 sc3 &>> log/pipeline-3.log;
        )&
    done
done

for scale_method in combat-gmm-2; do
    for cluster_method in sc3; do
        if [ ! -d results/cluster/${scale_method}-${cluster_method}/region ]; then
            mkdir results/cluster/${scale_method}-${cluster_method}/region
        fi

        if [ ! -d results/dimension_reduction/${scale_method}-${cluster_method} ]; then
            mkdir results/dimension_reduction/${scale_method}-${cluster_method}
        fi

        if [ ! -d draw_genes/${scale_method}-tsne ]; then
            mkdir draw_genes/${scale_method}-tsne
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
        wait $!

        # homer
        for region in regions; do
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

echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-3 finished."
