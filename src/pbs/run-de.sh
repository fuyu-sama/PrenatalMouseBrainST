PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

cd $HOME/workspace/mouse-brain-full/log
qsub << EOF
    #PBS -N DE
    #PBS -l nodes=1:ppn=1

    cd $HOME/workspace/mouse-brain-full

    Rscript src/R/run-de-region.R
    ${PYTHON_PATH} src/py/draw-de-region.py
    for region in cortex hypothalamus thalamus olfactory hippocampus
    do
        ${PYTHON_PATH} src/py/run-id-transfer.py \
            results/DE/region-specific/UP-\${region}.csv
    done
EOF
