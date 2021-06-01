for idx in P0B
do
    sleep 1
    cd $HOME/workspace/mouse-brain-full/spaceranger/log
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
    qsub << EOF
        #PBS -N ${idx}
        #PBS -l nodes=1:ppn=10
        #PBS -l walltime=240:00:00

        cd $HOME/workspace/mouse-brain-full/spaceranger
        source spaceranger-1.2.2.sh

        spaceranger count --id=${idx} \
                        --transcriptome=$HOME/refgenome/spaceranger/refdata-gex-mm10-2020-A \
                        --fastqs=$HOME/workspace/mouse-brain-full/Data/${idx} \
                        --image=$HOME/workspace/mouse-brain-full/Data/HE/${slide}-${idx}.tif \
                        --sample=${idx} \
                        --slide=${slide} \
                        --slidefile=$HOME/workspace/mouse-brain-full/spaceranger/slides/${slide}.gpr \
                        --area=${area} \
                        --localcores=10
        spaceranger mat2csv ${idx}/outs/filtered_feature_bc_matrix ${idx}/outs/filtered_feature_bc_matrix/${idx}.csv
EOF
done