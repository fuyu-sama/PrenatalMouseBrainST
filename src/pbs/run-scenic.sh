cd $HOME/workspace/mouse-brain-full
scale_method="combat-qq-logcpm"
cluster_method="sc3"
n_workers=20

source venv/bin/activate

if [ ! -d results/SCENIC/${scale_method}-${cluster_method} ]; then
    mkdir -p results/SCENIC/${scale_method}-${cluster_method}
fi

pyscenic grn \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}-${cluster_method}/adj.csv \
    results/SCENIC/${scale_method}-${cluster_method}/mean.csv \
    Data/SCENIC/mm_mgi_tfs.txt

pyscenic ctx \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}-${cluster_method}/regulons.csv \
    --annotations_fname Data/SCENIC/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
    --expression_mtx_fname results/SCENIC/${scale_method}-${cluster_method}/mean.csv  \
    results/SCENIC/${scale_method}-${cluster_method}/adj.csv \
    Data/SCENIC/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
    Data/SCENIC/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

pyscenic aucell \
    --num_workers ${n_workers} \
    -o results/SCENIC/${scale_method}-${cluster_method}/auc.csv \
    results/SCENIC/${scale_method}-${cluster_method}/mean.csv \
    results/SCENIC/${scale_method}-${cluster_method}/regulons.csv 
