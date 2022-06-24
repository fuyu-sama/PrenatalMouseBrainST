#PBS -N pipeline-3
#PBS -l nodes=1:ppn=40
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
    combat-Ai-0_500
    combat-Ai-500_1000
    combat-Ai-1000_1500
    combat-Ai-1500_2000
    combat-Ai-2000_2500
)

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
            Rscript src/R/run-rctd-${region}.R ${idx} &>> log/pipeline-3.log;
            ${PYTHON_PATH} src/py/draw-rctd.py \
                ${idx} ${region} combat-Ai-500union &>> log/pipeline-3.log;
        )&
    done
done

# %% dimension redunction
for scale_method in ${scale_methods[@]}; do
    if [ ! -d results/dimension_reduction/${scale_method}-sc3 ]; then
        mkdir results/dimension_reduction/${scale_method}-sc3
    fi
    (${PYTHON_PATH} src/py/run-dim_redu.py \
        ${scale_method} &>> log/pipeline-3.log)&
done

# %%
wait
echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-3 finished."
