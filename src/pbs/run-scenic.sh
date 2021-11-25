cd $HOME/workspace/mouse-brain-full
scale_method="combat-qq-logcpm"
n_workers=20

source venv/bin/activate

if [ ! -d results/SCENIC/${scale_method} ]; then
    mkdir -p results/SCENIC/${scale_method}
fi

pyscenic grn \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}/adj.csv \
    Data/scale_df/${scale_method}/full-${scale_method}.T.csv \
    Data/SCENIC/mm_mgi_tfs.txt

pyscenic ctx \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}/regulons.csv \
    --annotations_fname Data/SCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname Data/scale_df/${scale_method}/full-${scale_method}.T.csv  \
    results/SCENIC/${scale_method}/adj.csv \
    Data/SCENIC/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
    Data/SCENIC/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

pyscenic aucell \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}/auc.csv \
    Data/scale_df/${scale_method}/full-${scale_method}.T.csv \
    results/SCENIC/${scale_method}/regulons.csv 
