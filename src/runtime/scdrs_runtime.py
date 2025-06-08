#for running time analysis
#parameters:
# 1. anndata file path
# 2. gene set directory (all gene sets should contain .gs as the suffix)
# 3. result file path
import os
os.environ["OMP_NUM_THREADS"] = "1"

import pandas as pd
import scdrs
import argparse
import scanpy as sc
import time
import fnmatch
from joblib import dump, load

working_path = "seismic-analysis" #make it your working directory path

if __name__ == "__main__":
    os.chdir(working_path) #change it as your own directory

    #parse the argument
    #input parameters: 
    #1: expression file: h5ad
    #2: gene set directory
    #3: where to store the result file
    parser = argparse.ArgumentParser(description="Input parameters are: anndata, gene set directory, result file name")

    parser.add_argument('anndata_path', type=str, help='Path of the anndata file')
    parser.add_argument('gene_set_dir', type=str, help='Path of the gene set (.gs)')
    parser.add_argument('res_file_path', type=str, help='Path of the running time')

    args = parser.parse_args()

    #load expression file
    adata = scdrs.util.load_h5ad(args.anndata_path, flag_filter_data=True, flag_raw_count=True)

    files = os.listdir(args.gene_set_dir)
    gs_file_all = [file for file in files if fnmatch.fnmatch(file, '*.gs')]

    start = time.time()

    scdrs.preprocess(adata)
    sc.pp.pca(adata, n_comps=20) #default parameters
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20) #default paramters

    end = time.time()

    #processing time
    processing_time = end-start

    #calculating score time
    #group analysis time
    score_time = []
    group_time = []

    for gs_file in gs_file_all:
        
        gs_file = args.gene_set_dir + "/" + gs_file
        gs = scdrs.util.load_gs(gs_file, src_species="hsapiens", dst_species="hsapiens", to_intersect=adata.var_names)
        gene_list = [ value for key, value in gs.items()][0][0]
        gene_weight = [ value for key, value in gs.items()][0][1]
        
        time0 = time.time()
        df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000, return_ctrl_raw_score=False, return_ctrl_norm_score=True)
        time1 = time.time()
        group_res = scdrs.method.downstream_group_analysis(adata=adata, df_full_score= df_res, group_cols=["cell_ontology_class"])
        time2 = time.time()
        
        score_time.append(time1-time0)
        group_time.append(time2-time1)
    
    
    dump([processing_time, score_time, group_time], args.res_file_path)

