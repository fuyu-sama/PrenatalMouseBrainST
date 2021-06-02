PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full/SCT/log

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0B P0A1 P0A2
do
    sleep 1
    qsub << EOF
        #PBS -N ${idx}-SCTransform
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full/SCT

        source hdf5-1.12.0.sh
        Rscript $HOME/workspace/mouse-brain-full/src/R/run-SCTransform.R ${idx}
EOF
done
