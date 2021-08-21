#PBS -N pipeline-1
#PBS -l nodes=1:ppn=40
#PBS -l walltime=240:00:00

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2 P0B
do
    (
        cd $HOME/workspace/mouse-brain-full/Data/${idx};
        java -jar $Trimmomatic/trimmomatic-0.36.jar SE -threads 5 ${idx}_S1_L001_R1_001.fastq.gz ${idx}_S1_L001_R1_001.trimmed.fastq.gz CROP:28;
        rm ${idx}_S1_L001_R1_001.fastq.gz;
        mv ${idx}_S1_L001_R1_001.trimmed.fastq.gz ${idx}_S1_L001_R1_001.fastq.gz
    )&
done
wait

for idx in E135A E135B E155A E155B E165A E165B E175A1 E175A2 E175B P0A1 P0A2 P0B
do
    if [ $idx == "E135A" ]
    then
    area="A1"
    slide="V10M17-100"
    elif [ $idx == "E155A" ]
    then
    area="B1"
    slide="V10M17-100"
    elif [ $idx == "E165A" ]
    then
    area="C1"
    slide="V10M17-100"
    elif [ $idx == "P0B" ]
    then
    area="D1"
    slide="V10M17-100"
    elif [ $idx == "E175A1" ]
    then
    area="A1"
    slide="V10M17-101"
    elif [ $idx == "E175A2" ]
    then
    area="B1"
    slide="V10M17-101"
    elif [ $idx == "P0A1" ]
    then
    area="C1"
    slide="V10M17-101"
    elif [ $idx == "P0A2" ]
    then
    area="D1"
    slide="V10M17-101"
    elif [ $idx == "E135B" ]
    then
    area="A1"
    slide="V10M17-085"
    elif [ $idx == "E155B" ]
    then
    area="B1"
    slide="V10M17-085"
    elif [ $idx == "E165B" ]
    then
    area="C1"
    slide="V10M17-085"
    elif [ $idx == "E175B" ]
    then
    area="D1"
    slide="V10M17-085"
    fi

    cd $HOME/workspace/mouse-brain-full/spaceranger
    (
        source spaceranger-1.2.2.sh;
        spaceranger count \
            --id=${idx} \
            --transcriptome=$HOME/refgenome/spaceranger/refdata-gex-mm10-2020-A \
            --fastqs=$HOME/workspace/mouse-brain-full/Data/${idx} \
            --image=$HOME/workspace/mouse-brain-full/Data/HE/${slide}-${idx}.tif \
            --sample=${idx} \
            --slide=${slide} \
            --slidefile=$HOME/workspace/mouse-brain-full/spaceranger/slides/${slide}.gpr \
            --area=${area} \
            --localcores=10 >>log/pipeline-1.log;
        spaceranger mat2csv \
            ${idx}/outs/filtered_feature_bc_matrix \
            ${idx}/outs/filtered_feature_bc_matrix/${idx}.csv
    )&
done
wait

echo "[`date +%Y.%m.%d\ %H:%M:%S`] Pipeline-1 finished."
