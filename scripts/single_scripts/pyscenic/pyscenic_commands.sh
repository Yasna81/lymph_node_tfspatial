#!/bin/bash
#-----------------------------------------------
#pyscenic workflow summary
#----------------------------------------------
#adjust the pathes-file names before running!
﻿#first activate your env and install the pyscenic 
#conda create -n scenic-set python=3.10 -y
conda activate scenic-set
conda install -c conda-forge -c bioconda pyscenic -y 
# make sure you have installed every thing properly 
pyscenic grn --version
pyscenic –help
#use the transposed version of expression matrix by the python script (transpos.py) > cells in rows and genes in cols
#step 1 : expr_matrix + tf list for human > adjacencies (co expression links between TFs and target genes)
pyscenic grn expr_matrix_transposed.csv allTFs_hg38.txt -o adjacencies.tsv --num_workers 4
#step 2 : adjacencies + motif database > rgc.csv (which TF-target links have shared DNA motif/biologically whcih tf binds to dna/)
pyscenic ctx adjacencies.tsv \hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \--annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \--expression_mtx_fname expr_matrix_transposed.csv \--mode "dask_multiprocessing" \--output reg.csv \--num_workers 1 \--top_n_targets 200 \--top_n_regulators 500 \--mask_dropouts \--mode custom_m
#step 3 :regulon activity per cell, which TF which cell > auc_mtx.cv
pyscenic aucell expr_matrix_transposed.csv reg.csv -o auc_mtx.csv --num_workers 1
#step 4 : followed by umap and violin plots of clusters/tfs in R+suerat object :)
