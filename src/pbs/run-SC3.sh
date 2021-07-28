PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

scale_method=scMerge
for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0B P0A1 P0A2
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/log
    qsub << EOF
        #PBS -N ${scale_method}-SC3-${idx}
        #PBS -l nodes=comput3:ppn=5

        cd $HOME/workspace/mouse-brain-full

        Rscript $HOME/workspace/mouse-brain-full/src/R/run-SC3.R ${idx} ${scale_method}
        ${PYTHON_PATH} $HOME/workspace/mouse-brain-full/src/py/draw-SC3.py ${idx} ${scale_method}

EOF
done
