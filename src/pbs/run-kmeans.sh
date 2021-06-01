PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0B P0A1 P0A2 
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/logcpm/log
    qsub << EOF
        #PBS -N ${idx}-kmeans
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full/logcpm

        ${PYTHON_PATH} $HOME/workspace/mouse-brain-full/src/py/run-kmeans.py ${idx}

EOF
done