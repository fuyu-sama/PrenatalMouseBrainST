WORKDIR=$HOME/workspace/mouse-brain-full
cd $WORKDIR

ulimit -n 10240
ulimit -v 33170449147
NUMBA_CACHE_DIR=$WORKDIR/log

dataDir=$WORKDIR/Data/Stereo-seq
fq1=$dataDir/DP8450004027TR_L01_57_1.fq.gz
fq2=$dataDir/DP8450004027TR_L01_57_2.fq.gz
gtf=$HOME/refgenome/SAW/spaceranger-refdata-gex-mm10-2020-A/genes/genes.gtf
maskfile=$dataDir/SS200000977TL_A1.barcodeToPos.h5
refindex=$HOME/refgenome/SAW/spaceranger-refdata-gex-mm10-2020-A/STAR_SJ100
imgrecord=$dataDir/QCImgUpload/SS200000977TL_A1_20221028_140744_1.0.8.json
imgcompress=$dataDir/QCImgUpload/SS200000977TL_A1_20221028_140744_1.0.8.tar.gz
imgdir=$dataDir/QCImgUpload
for n in `seq 58 1 64`; do
    fq1=${fq1},${dataDir}/DP8450004027TR_L01_${n}_1.fq.gz
    fq2=${fq2},${dataDir}/DP8450004027TR_L01_${n}_2.fq.gz
done
export SINGULARITY_BIND=$dataDir,$outDir

if false; then
    outDir=$WORKDIR/SAW.5.1.3.with.img
    if [[ ! -d ${outDir} ]]; then
        mkdir $outDir
    fi
    bash $HOME/workspace/SAW/SAW/script/stereoPipeline.sh \
        -genomeSize 25 \
        -splitCount 1 \
        -maskFile $maskfile \
        -fq1 $fq1 \
        -fq2 $fq2 \
        -speciesName mouse \
        -tissueType brain \
        -refIndex $refindex \
        -annotationFile $gtf \
        -sif $HOME/workspace/SAW/SAW_v5.1.3.sif \
        -imageRecordFile $imgrecord\
        -imageCompressedFile $imgcompress \
        -threads 10 \
        -doCellBin Y \
        -outDir $outDir
fi

if true; then
    outDir=$WORKDIR/SAW.4.1.0.with.img
    if [[ ! -d ${outDir} ]]; then
        mkdir $outDir
    fi
    bash $HOME/workspace/SAW/SAW/script/stereoRun_multiLane_v4.1.0.sh \
        -m $maskfile \
        -1 $fq1 \
        -2 $fq2 \
        -g $refindex \
        -a $gtf \
        -o $outDir \
        -i $imgdir \
        -t 10 \
        -s $HOME/workspace/SAW/SAW_v4.1.0.sif \
        -c 25
fi
