#file to do null simulation
import os
import pandas as pd
import scdrs
import argparse
import scanpy as sc
from functools import partial
from joblib import dump, load
from multiprocessing import Pool


def ass_mcz(rand_cell_idx, score_df, anndata):
    anndata.obs["cell_ontology_class"] = "others"
    anndata.obs.iloc[rand_cell_idx,anndata.obs.columns.get_loc("cell_ontology_class")] = "fake_cell_type"
    group_res = scdrs.method.downstream_group_analysis(adata=anndata, df_full_score= score_df, group_cols=["cell_ontology_class"])
    return group_res["cell_ontology_class"].iloc[0]["assoc_mcz"]


def comp_seed_table_ass_idx(seed_table, score_df, anndata, idx):
    if idx % 10 == 0:
        print("Processing index:"+str(idx))
    return ass_mcz(seed_table.iloc[idx]-1, score_df, anndata) #1 based but not 0 based

def comp_seed_table_ass(seed_table, score_df, anndata, num_core):
    f = partial(comp_seed_table_ass_idx, seed_table, score_df,anndata)
    with Pool(num_core) as p:
        res = p.map(f, range(0,10000))
    return res

if __name__ == "__main__":
    os.chdir("/grain/ql29/seismic-analysis") #set the working directory
 
    #parse the argument
    #input parameters: 
    #1: anndata loading
    #2: seed table
    #3: gene set file
    #4 num of cores
    #5 output files
    parser = argparse.ArgumentParser(description="Input parameters are: anndata, seed table, z-score file, number of cores and output files")

    parser.add_argument('anndata_path', type=str, help='Path of the anndata file')
    parser.add_argument('seed_table_path', type=str, help='Path of the seed file: for sampling')
    parser.add_argument('gs_path', type=str, help='Path of the gene set file')
    parser.add_argument('num_cores', type=int, help='Number of cores')
    parser.add_argument('output_path', type=str, help='Path for the output result file')

    args = parser.parse_args()

    #load seed file
    seed_table = pd.read_csv(args.seed_table_path, sep="\t")

    #load data and process the data
    adata = scdrs.util.load_h5ad(args.anndata_path, flag_filter_data=True, flag_raw_count=True)
    gs = scdrs.util.load_gs(args.gs_path, src_species="hsapiens", dst_species="mmusculus", to_intersect=adata.var_names)
    scdrs.preprocess(adata)
    gene_list = [ value for key, value in gs.items()][0][0]
    gene_weight = [ value for key, value in gs.items()][0][1]
    sc.pp.pca(adata, n_comps=20) #default parameters
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20) #default paramters

    #compute score
    df_res = scdrs.score_cell(adata, gene_list, gene_weight=gene_weight, n_ctrl=1000, return_ctrl_raw_score=False, return_ctrl_norm_score=True)
    ass_mcz_res_all = comp_seed_table_ass(seed_table, df_res, adata, args.num_cores)
    dump(ass_mcz_res_all,args.output_path[:-3]+"joblib") #dump temp object
    out_df = pd.DataFrame({"index": range(1,10001), "assoc_mcz":ass_mcz_res_all})
    out_df.to_csv(args.output_path, sep="\t",index=False)
