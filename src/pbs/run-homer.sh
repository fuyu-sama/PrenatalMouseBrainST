PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

scale_method=combat
for region in cortex hypothalamus thalamus olfactory hippocampus
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/log
    qsub << EOF
        #PBS -N ${region}-${scale_method}-homer
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full

        source homer-4.11.sh
        findMotifsGenome.pl \
            results/DE/${scale_method}/region-specific/motifInputs/${region}.bed \
            mm10 results/motifResults/${scale_method}/${region}
EOF
done
