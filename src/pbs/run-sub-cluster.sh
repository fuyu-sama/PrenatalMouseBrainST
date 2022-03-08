PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
    if [ ! -d results/5/${idx} ]; then
        mkdir results/5/${idx}
    fi
    ${PYTHON_PATH} src/test-scripts/run-gene-cluster.py cpm ${idx}
done
