PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

for region in cortex hippocampus hypothalamus thalamus
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/logcpm/log
    qsub << EOF
        #PBS -N ${region}-homer
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full/logcpm/DE/region-specific

        source homer-4.11.sh
        ${PYTHON_PATH} $HOME/workspace/mouse-brain-full/src/py/run-id-transfer.py UP-12-${region}.csv
        findMotifs.pl UP-12-${region}.csv.out mouse motifResults/${region} \
            -start -1000 -end 1000
EOF
done
