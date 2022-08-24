#PBS -N sc3-all
#PBS -l nodes=comput2:ppn=20
#PBS -l walltime=2400:00:00

# %% environment config
PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python
idx_full=(
    E135A E135B 
    E155A E155B 
    E165A E165B 
    E175A1 E175A2
    P0A1 P0A2
)

cd $HOME/workspace/mouse-brain-full

# %%
scale_method="combat"
if true; then
    for idx in ${idx_full[@]}; do
        gene_list=Data/gene-lists/${idx}-logcpm-0.99.csv
        suffix=${idx}-logcpm-0.99
        if [ ! -d results/cluster/${scale_method}-${suffix}-sc3 ]; then
            mkdir results/cluster/${scale_method}-${suffix}-sc3
            mkdir results/cluster/${scale_method}-${suffix}-sc3/pattern
            mkdir results/cluster/${scale_method}-${suffix}-sc3/separate
            mkdir results/cluster/${scale_method}-${suffix}-sc3/together
        fi
        (
            Rscript src/R/run-sc3.R \
                ${idx} ${scale_method} ${gene_list} ${suffix} &>> log/sc3.log;
            ${PYTHON_PATH} src/py/draw-sc3.py \
                ${idx} ${scale_method} ${suffix} &>> log/sc3.log;
        )&
    done
    wait
fi

# %%
step_len=500
end_point=`expr 2000 + $step_len`
scale_method=combat
if false; then
    for n in `seq 0 ${step_len} ${end_point}`; do
        gene_list=logcpm-0.99-Ai-${n}_`expr ${n} + ${step_len}`
        gene_list=logcpm-hotspot-6-logcpm-0.99-Ai-${n}_`expr ${n} + ${step_len}`
        if [ ! -d results/cluster/${scale_method}-${gene_list}-sc3 ]; then
            mkdir results/cluster/${scale_method}-${gene_list}-sc3
            mkdir results/cluster/${scale_method}-${gene_list}-sc3/pattern
            mkdir results/cluster/${scale_method}-${gene_list}-sc3/separate
            mkdir results/cluster/${scale_method}-${gene_list}-sc3/together
        fi
        for idx in ${idx_full[@]}; do
            (
                Rscript src/R/run-sc3.R \
                    ${idx} ${scale_method} \
                    Data/gene-lists/${idx}-${gene_list}.csv \
                    ${gene_list} &>> log/sc3.log;
                ${PYTHON_PATH} src/py/draw-sc3.py \
                    ${idx} ${scale_method} ${gene_list} &>> log/sc3.log;
            )&
        done
        wait
        (${PYTHON_PATH} src/py/run-dim_redu.py ${gene_list})&
    done
    wait
fi
