#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path


sample_id = sys.argv[1]
output_file = sys.argv[2]
input_files = sys.argv[3:]

# Load and concatenate the CSVs by columns
dfs = [pd.read_csv(f, dtype=str) for f in input_files]
combined_df = pd.concat(dfs, axis=1)

# Insert Sample_ID as first column
combined_df.insert(0, "Sample_ID", sample_id)

# Reorder columns: Sample_ID → all QC columns → the rest of the columns
fixed_first = ['Sample_ID']
qc_cols = [
    col for col in combined_df.columns
    if col != 'Sample_ID' and combined_df[col].dropna().isin(['PASS', 'FAIL']).all()
]
other_cols = [col for col in combined_df.columns if col not in fixed_first + qc_cols]
combined_df = combined_df[fixed_first + qc_cols + other_cols]

# Save to output
combined_df.to_csv(output_file, index=False)