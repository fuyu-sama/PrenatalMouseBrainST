for idx in P0B
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/Data/${idx}
    qsub << EOF
        #PBS -N ${idx}-trim
        #PBS -l nodes=comput3:ppn=5
        #PBS -l walltime=24:00:00

        cd $HOME/workspace/mouse-brain-full/Data/${idx}

        java -jar $Trimmomatic/trimmomatic-0.36.jar SE -threads 5 ${idx}_S1_L001_R1_001.fastq.gz ${idx}_S1_L001_R1_001.trimmed.fastq.gz CROP:28
        rm ${idx}_S1_L001_R1_001.fastq.gz
        mv ${idx}_S1_L001_R1_001.trimmed.fastq.gz ${idx}_S1_L001_R1_001.fastq.gz
EOF
done