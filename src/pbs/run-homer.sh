PYTHON_PATH=$HOME/workspace/mouse-brain-full/venv/bin/python

for region in cortex hypothalamus thalamus olfactory hippocampus
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/log
    qsub << EOF
        #PBS -N ${region}-homer
        #PBS -l nodes=1:ppn=1

        cd $HOME/workspace/mouse-brain-full

        source homer-4.11.sh
        findMotifsGenome.pl results/DE/region-specific/motifInputs/${region}.bed \
            mm10 results/motifResults/${region}
EOF
done
