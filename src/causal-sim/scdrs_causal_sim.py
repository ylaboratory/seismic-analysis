# get causal simulation result for each cell type
# input arguments for the script are:
# 1. parameter data frame path
# 2. output directory
# 3. output header column
# 4. target cell type pattern (regular expression to extract the target cell type)
# 5. extract target cell (extract the statistics of the target cell type or not.
#    Any value that does not match any cell types will return 0.)
# 6. number of cores to use

import os
import argparse
import multiprocessing as mp
from multiprocessing import Pool
from functools import partial

import pandas as pd
import psutil
import scdrs
import scanpy as sc
import anndata as ad
import numpy as np
import threadpoolctl  # control PCA thread
from joblib import dump

# set environment variables
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

working_path="seismic-analysis" #replace to your own working directory

def worker_init():
    thread_id = int(mp.current_process().name.split('-')[1])

    # set environment variables
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

    # control the number of threads for PCA
    threadpoolctl.threadpool_limits(limits=1)

    # thread id to cpu id
    thread_id = thread_id % 72
    cpu_id = [thread_id * 2, thread_id * 2 + 1]
    p = psutil.Process()
    p.cpu_affinity(cpu_id)


# adapt from scdrs method for customized anndata preparation
def prepare_adata(adata_obj, replace_mat, new_cell_anno, group="cell_ontology_class"):
    # label cell type where the cellid matches the colnames of the target cell type
    perturbed_cell_id = replace_mat.columns
    renamed_cell_id = new_cell_anno.index

    # raise error if any cellid is not in original anndata
    missing_cell_ids = (set(perturbed_cell_id) - set(adata_obj.obs["cellid"])).union(
        set(renamed_cell_id) - set(perturbed_cell_id)
    )
    if missing_cell_ids:
        raise ValueError(f"Cell IDs {missing_cell_ids} are not present in the original anndata object.")

    # get the unique cell type from new cell anno, which is a series of cell type with index equals to cell names
    new_cell_type = [i for i in new_cell_anno.unique() if i not in adata_obj.obs[group].cat.categories]
    perturbed_adata = adata_obj.copy()
    perturbed_adata.obs[group] = perturbed_adata.obs[group].cat.add_categories(new_cell_type)

    # subsititute the cell type where the cellid matches the colnames of the target cell type to new cell anno
    perturbed_adata.obs = perturbed_adata.obs.set_index("cellid")
    perturbed_adata.obs.loc[renamed_cell_id, group] = new_cell_anno.loc[renamed_cell_id]

    # add gene perturbation infromation
    perturbed_genes = replace_mat.index.values

    # matrix replacement
    row_indicies = [np.where(perturbed_adata.obs.index == i)[0][0] for i in perturbed_cell_id]
    col_indicies = pd.Index(perturbed_adata.var["symbol"]).get_indexer(perturbed_genes)
    perturbed_adata[row_indicies, col_indicies].X = np.array(replace_mat).T

    # filter out cell types with fewer than 20 cells
    cell_counts = perturbed_adata.obs[group].value_counts()
    valid_cell_types = cell_counts[cell_counts >= 20].index
    perturbed_adata = perturbed_adata[perturbed_adata.obs[group].isin(valid_cell_types)].copy()

    # filter cells
    sc.pp.filter_cells(perturbed_adata, min_genes=250)  # default parameter
    sc.pp.filter_genes(perturbed_adata, min_cells=50)  # default parameter

    # log transform
    sc.pp.normalize_per_cell(perturbed_adata, counts_per_cell_after=1e4)
    sc.pp.log1p(perturbed_adata)

    return perturbed_adata


def preprocess_scdrs(adata_obj):
    # control the number of threads for PCA
    threadpoolctl.threadpool_limits(1)

    # preprocessda
    scdrs.preprocess(adata_obj)
    sc.pp.pca(adata_obj, n_comps=20)  # default parameters
    sc.pp.neighbors(adata_obj, n_neighbors=15, n_pcs=20)  # default paramters
    return adata_obj


def get_causal_sim_res(adata_obj, gene_set, group="cell_ontology_class"):
    # gene list and weight
    gene_list = [value for key, value in gene_set.items()][0][0]
    gene_weight = [value for key, value in gene_set.items()][0][1]

    # score cells
    df_res = scdrs.score_cell(
        adata_obj, gene_list, gene_weight, n_ctrl=1000,
        return_ctrl_raw_score=False, return_ctrl_norm_score=True
    )
    group_res = scdrs.method.downstream_group_analysis(adata=adata_obj, df_full_score=df_res, group_cols=[group])

    return group_res[group]


def scdrs_file_to_res(
    adata_file, replace_mat_file, cell_anno_file, gs_file,
    output_path, target_cell_type_pattern, extract_target_cell,
    group_to_replace="new_cell_ontology_class", group="cell_ontology_class"
):
    # if the output file exist, then the results equals to the output file
    if os.path.exists(output_path):
        res = pd.read_csv(output_path, sep="\t", index_col=0)
        # filter index with target cell type pattern
        res = res[res.index.str.contains(target_cell_type_pattern)]

    else:
        threadpoolctl.threadpool_limits(1)
        # load data
        adata = ad.read(adata_file)
        replace_mat = pd.read_csv(replace_mat_file, sep=" ")
        cell_anno = pd.read_csv(cell_anno_file, sep="\t")
        gene_set = scdrs.util.load_gs(
            gs_file, src_species="hsapiens", dst_species="mmusculus", to_intersect=adata.var_names
        )

        # prepare cell annotation
        cell_anno = cell_anno.set_index("cellid")

        # get the new cell annotation: where index matches the column names of replace_mat
        new_cell_anno = cell_anno.loc[replace_mat.columns, group_to_replace]

        # prepare adata
        perturbed_adata = prepare_adata(adata, replace_mat, new_cell_anno, group)

        # preprocess adata
        perturbed_adata = preprocess_scdrs(perturbed_adata)

        # get causal simulation result
        res = get_causal_sim_res(perturbed_adata, gene_set, group)
        res.to_csv(output_path, sep="\t", index=True)

    # filter index with target cell type pattern
    res = res[res.index.str.contains(target_cell_type_pattern)]

    if res.shape[0] != 0 and extract_target_cell:
        return (list(res.index), list(res["assoc_mcz"]))
    else:
        return (["Error_ct"], [0])


def parameter_df_to_res(parameter_df, output_header_col, target_cell_type_pattern, extract_target_cell, idx):
    # use arguments plus an index (indicating the row of the parameter_df) to get the results of current case
    adata_file = parameter_df.loc[idx, "expr_h5ad_file"]
    replace_mat_file = parameter_df.loc[idx, "perturbed_mat_file"]
    gs_file = parameter_df.loc[idx, "gs_scdrs_file"]
    cell_anno_file = parameter_df.loc[idx, "cell_anno_file"]
    output_file_path = parameter_df.loc[idx, output_header_col] + ".scdrs_res.csv"

    res = scdrs_file_to_res(
        adata_file, replace_mat_file, cell_anno_file, gs_file, output_file_path,
        target_cell_type_pattern=target_cell_type_pattern, extract_target_cell=extract_target_cell
    )

    if idx % 10 == 0:
        print("Processing index:" + str(idx))
        print("The file is already in " + output_file_path)

    return res


def main():
    """Main function to run the script."""
    os.chdir(working_path)  # set the working directory

    parser = argparse.ArgumentParser(
        description=("Input parameters are: parameter data frame, output directory, "
                     "organized_by column and number of cores")
    )
    
    parser.add_argument('parameter_df_path', type=str, help='Path of the parameter data frame')
    parser.add_argument('final_res_path', type=str, help='Path of the output file')
    parser.add_argument(
        'output_header_col', type=str,
        help='In the output directory, organize the output by the directories named after this column'
    )
    parser.add_argument('target_ct_pattern', type=str, help='Regular expression to extract the target cell type')
    parser.add_argument(
        'extract_target_cell', type=str,
        help=("Extract the statistics of the target cell type or not. "
              "Any value that does not match any cell types will return 0.")
    )
    parser.add_argument('num_cores', type=int, help='Number of cores')

    args = parser.parse_args()

    # read in df
    parameter_df = pd.read_csv(args.parameter_df_path, sep=" ")

    output_dir_all = {os.path.dirname(i) for i in parameter_df[args.output_header_col].unique()}
    for value in output_dir_all:
        if not os.path.exists(value):
            os.makedirs(value, exist_ok=True)

    f = partial(
        parameter_df_to_res, parameter_df, args.output_header_col,
        args.target_ct_pattern, bool(args.extract_target_cell)
    )

    with Pool(args.num_cores, initializer=worker_init) as p:
        asso_mcz_all = p.map(f, range(0, parameter_df.shape[0]))

    dump(asso_mcz_all, args.final_res_path + "/scdrs_all_res.joblib")

    # create a data frame using asso_mcz_all, with column name
    if args.extract_target_cell and len(asso_mcz_all) > 0 and asso_mcz_all[0][0] != "Error_ct":
        all_ct = set()
        for cur_res in asso_mcz_all:
            all_ct.update(cur_res[0])
        if all_ct == set(asso_mcz_all[0][0]):
            out_df = pd.DataFrame({"index": range(1, parameter_df.shape[0] + 1)})
            for i, ct in enumerate(asso_mcz_all[0][0]):
                out_df[ct + "_assoc_mcz"] = [res[1][i] if i < len(res[1]) else None for res in asso_mcz_all]
            out_df.to_csv(os.path.join(args.final_res_path, "scdrs_all_res"), index=False)


if __name__ == "__main__":
    main()
