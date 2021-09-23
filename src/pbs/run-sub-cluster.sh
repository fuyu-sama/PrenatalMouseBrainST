#PBS -N sub-sc3
#PBS -l nodes=1:ppn=20
#PBS -l walltime=240:00:00

PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
    (
    if [ ! -d results/cluster/combat-sc3/pattern/sub-clusters ]; then
        mkdir results/cluster/combat-sc3/pattern/sub-clusters
    fi

    if [ ! -d results/cluster/combat-sc3/together/sub-clusters ]; then
        mkdir results/cluster/combat-sc3/together/sub-clusters
    fi

    Rscript src/R/run-sub-sc3.R ${idx} combat
    ${PYTHON_PATH} src/py/draw-sub-sc3.py ${idx} combat
    )&
done
wait
