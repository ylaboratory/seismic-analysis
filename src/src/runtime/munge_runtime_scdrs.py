#import libraries
import os
import pandas as pd
from joblib import load

os.chdir("/grain/ql29/seismic-analysis")

#load all files 
sample_size = ["10k", "25k", "50k", "100k", "150k", "200k", "250k", "300k"]
ds = ["ds_1", "ds_2", "ds_3", "ds_4", "ds_5"]
all_df = []

for s in sample_size:
    for d in ds:
        runtime_case = load(f"results/runtime/{s}/{d}/{d}.scdrs.joblib")
        df = pd.DataFrame({"processing_time": runtime_case[0],
              "score_time": runtime_case[1],
              "association_time": runtime_case[2],
              "size":s,
              "ds":d})
        all_df.append(df)

#concatenate all dataframes
all_df = pd.concat(all_df, ignore_index=True)

#print out as a cvs
all_df.to_csv("results/runtime/scDRS/runtime.txt", index=False, sep="\t")