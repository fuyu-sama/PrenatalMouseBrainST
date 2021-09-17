#PBS -N RCTD
#PBS -l nodes=1:ppn=20
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

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
        Rscript src/R/run-rctd-${region}.R ${idx} combat sc3;

        exit_status=$?
        if [ $exit_status -eq 0 ]; then
            ${PYTHON_PATH} src/py/draw-rctd.py ${idx} ${region};
        fi
        )&
    done
done
wait
