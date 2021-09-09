#PBS -N pipeline-2
#PBS -l nodes=1:ppn=10
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
    (
        if [ ! -d results/RCTD/${idx} ]; then
            mkdir results/RCTD/${idx}
            mkdir results/RCTD/${idx}/results
            mkdir results/RCTD/${idx}/weights
        fi

        source hdf5-1.12.0.sh
        Rscript src/R/run-rctd.R ${idx} combat sc3;
        ${PYTHON_PATH} src/py/draw-rctd.py ${idx};
    )&
done
wait
