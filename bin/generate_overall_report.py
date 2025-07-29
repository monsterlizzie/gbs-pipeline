#! /usr/bin/env python3

# Generate overall report based on sample reports and columns specified by COLUMNS_BY_CATEGORY and ARIBA metadata

import sys
from itertools import chain
import pandas as pd
import glob


# Columns to include in output, organized by category
COLUMNS_BY_CATEGORY = {
    'IDENTIFICATION': ['Sample_ID'],
    'QC': ['Read_QC', 'Assembly_QC', 'Mapping_QC', 'Taxonomy_QC', 'Overall_QC'],
    'READ': ['Bases'],
    'ASSEMBLY': ['Contigs#', 'Assembly_Length', 'Seq_Depth'],
    'MAPPING': ['Ref_Cov_%', 'Het-SNP#'],
    'TAXONOMY': ['S.agalactiae_%', 'Top_Non-Agalactiae_Species', 'Top_Non-agalactiae_Species_%']
}

OUTPUT_COLUMNS = list(chain.from_iterable(COLUMNS_BY_CATEGORY.values()))


def get_df_output(input_pattern, output_columns):
    df_manifest = pd.DataFrame(columns=output_columns)

    dfs = [df_manifest]
    reports = glob.glob(input_pattern)
    for report in reports:
        df = pd.read_csv(report, dtype=str)  
        df.replace("", pd.NA, inplace=True)  
        dfs.append(df)

    df_output = pd.concat(dfs, ignore_index=True).sort_values(by='Sample_ID')

    # Ensure all expected columns exist
    for col in output_columns:
        if col not in df_output.columns:
            df_output[col] = pd.NA

    df_output = df_output[output_columns]

    # Fill missing fields in PASS samples with "MODULE FAILURE"
    df_output.loc[df_output["Overall_QC"] == "PASS"] = (
        df_output.loc[df_output["Overall_QC"] == "PASS"].fillna(value="MODULE FAILURE")
    )

    # Fill all other missing values with "NA"
    df_output.fillna("NA", inplace=True)

    return df_output


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: generate_overall_qc_report.py <input_pattern> <output_file>")

    input_pattern = sys.argv[1]
    output_file = sys.argv[2]

    df_output = get_df_output(input_pattern, OUTPUT_COLUMNS)
    df_output.to_csv(output_file, index=False)


if __name__ == "__main__":
    main()