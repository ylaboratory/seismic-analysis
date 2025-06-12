# import libraries
import os
import pandas as pd
from joblib import load

working_path = "seismic-analysis"  # make it your working directory path

os.chdir(working_path)

# A list to hold each DataFrame from the loop
list_of_dfs = []

# load all files
sample_size = ["10k", "25k", "50k", "100k", "150k", "200k", "250k", "300k", "400k", "500k"]
ds = ["ds_1", "ds_2", "ds_3", "ds_4", "ds_5"]

for s in sample_size:
    for d in ds:
        runtime_case = load(f"results/runtime/{s}/{d}/{d}.scdrs.joblib")
        df = pd.DataFrame({"processing_time": runtime_case[0],
                           "score_time": runtime_case[1],
                           "association_time": runtime_case[2],
                           "size": s,
                           "ds": d})
        list_of_dfs.append(df)

# concatenate all dataframes into a new variable
final_df = pd.concat(list_of_dfs, ignore_index=True)

# print out as a csv
final_df.to_csv("results/runtime/scDRS/runtime.txt", index=False, sep="\t")
