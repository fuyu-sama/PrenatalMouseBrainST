PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2; do
    for region in hypothalamus cortex thalamus; do
            (${PYTHON_PATH} src/py/run-sub-cluster.py combat-zjq sc3 ${idx} ${region})&
    done
    wait
done