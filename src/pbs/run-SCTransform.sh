PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full/log

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0B P0A1 P0A2
do
    sleep 1
    qsub << EOF
        #PBS -N SCTransform-${idx}
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full

        Rscript $HOME/workspace/mouse-brain-full/src/R/run-SCTransform.R ${idx}
EOF
done
