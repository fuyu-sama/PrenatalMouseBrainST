#PBS -N pipeline-3
#PBS -l nodes=1:ppn=40
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
    logcpm
    cpm
    raw
)
regions=(cortex hippocampus hypothalamus thalamus region_a)

cd $HOME/workspace/mouse-brain-full


# %% RCTD
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
                ${idx} combat-logcpm-1000 sc3 &>> log/pipeline-3.log;

            ${PYTHON_PATH} src/py/draw-rctd.py \
                ${idx} ${region} combat-logcpm-1000 sc3 &>> log/pipeline-3.log;
            #${PYTHON_PATH} src/py/draw-rctd-IE.py \
                #${idx} ${region} combat-gmm-2 sc3 &>> log/pipeline-3.log;
        )&
    done
done

# %% draw regions
for scale_method in combat; do
    for cluster_method in sc3; do
        if [ ! -d results/cluster/${scale_method}-${cluster_method}/region ]; then
            mkdir results/cluster/${scale_method}-${cluster_method}/region
        fi
        (${PYTHON_PATH} src/py/draw-region.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done

# %% dimension redunction
for scale_method in combat; do
    for cluster_method in sc3; do
        if [ ! -d results/dimension_reduction/${scale_method}-${cluster_method} ]; then
            mkdir results/dimension_reduction/${scale_method}-${cluster_method}
        fi
        if [ ! -d draw_genes/${scale_method}-tsne ]; then
            mkdir draw_genes/${scale_method}-tsne
        fi
        (${PYTHON_PATH} src/py/run-dim_redu.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done
        
# %% draw hierarchical tree
for scale_method in combat; do
    for cluster_method in sc3; do
        (Rscript src/R/run-hierarchical.R \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done

# %% DE
for scale_method in combat; do
    for cluster_method in sc3; do
        if [ ! -d results/DE/${scale_method}-${cluster_method} ]; then
            mkdir results/DE/${scale_method}-${cluster_method}
            mkdir results/DE/${scale_method}-${cluster_method}/region-specific
            mkdir results/DE/${scale_method}-${cluster_method}/timepoint-specific
        fi

        (
            Rscript src/R/run-de-region.R \
                ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
            ${PYTHON_PATH} src/py/draw-de-region.py \
                ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
            Rscript src/R/run-de-timepoint.R \
                ${scale_method} ${cluster_method} &>> log/pipeline-3.log;
            ${PYTHON_PATH} src/py/draw-de-timepoint.py \
                ${scale_method} ${cluster_method} &>> log/pipeline-3.log
        )& de_pids+=($!)
    done
done

# %% homer
for scale_method in combat; do
    for cluster_method in sc3; do
        for region in ${regions[@]}; do
            wait "${de_pids[@]}"
            (
                source homer-4.11.sh;
                ${PYTHON_PATH} src/py/run-id-transfer.py \
                    results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv \
                    &>> log/pipeline-3.log;
                findMotifs.pl \
                    results/DE/${scale_method}-${cluster_method}/region-specific/UP-${region}.csv.out \
                    mouse \
                    results/motifResults/${scale_method}-${cluster_method}/${region} \
                    -start -1000 -end 1000 &>> log/pipeline-3.log
            )& homer_pids+=($!)
        done
        wait "${homer_pids[@]}"
        (${PYTHON_PATH} src/py/draw-homer.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
        (${PYTHON_PATH} src/py/draw-go.py \
            ${scale_method} ${cluster_method} &>> log/pipeline-3.log)&
    done
done

# %% gene cluster homer
n_gene_clusters=12
for idx in ${idx_full[@]}; do
    for i in {1..$n_gene_clusters}; do
        (
            source homer-4.11.sh;
            ${PYTHON_PATH} src/py/run-id-transfer.py \
                results/gene-cluster/${idx}/tables/${idx}-genes-${i}.csv \
                &>> log/pipeline-3.log
            findMotifs.pl \
                results/gene-cluster/${idx}/tables/${idx}-genes-${i}.csv.out \
                mouse \
                results/motifResults/gene-cluster/${idx}-${i} \
                -start -1000 -end 1000 &>> log/pipeline-3.log
        )& gene_cluster_homer_pids+=($!)
    done
done

wait "${gene_cluster_homer_pids[@]}"
for idx in ${idx_full[@]}; do
    (${PYTHON_PATH} src/py/draw-genes-homer.py \
        ${idx} ${n_gene_clusters} &>> log/pipeline-3.log)&
    (${PYTHON_PATH} src/py/draw-genes-go.py \
        ${idx} ${n_gene_clusters} &>> log/pipeline-3.log)&
done

wait
echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-3 finished."
