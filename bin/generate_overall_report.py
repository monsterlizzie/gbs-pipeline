#!/usr/bin/env python3

import sys, os, re, glob
import pandas as pd
from itertools import chain
from typing import List

# ---------- Columns in include ----------
COLUMNS_BY_CATEGORY = {
    'IDENTIFICATION': ['Sample_ID'],
    'QC': ['Read_QC', 'Assembly_QC', 'Mapping_QC', 'Taxonomy_QC', 'Overall_QC'],
    'READ': ['Bases'],
    'ASSEMBLY': ['Contigs#', 'Assembly_Length', 'Seq_Depth'],
    'MAPPING': ['Ref_Cov_%', 'Het-SNP#'],
    'TAXONOMY': ['S.agalactiae_%', 'Top_Non-agalactiae_Species', 'Top_Non-agalactiae_Species_%'],
}
QC_FIXED = list(chain.from_iterable(COLUMNS_BY_CATEGORY.values()))

VALID_ID = re.compile(r'^[A-Za-z0-9_.:-]+$')

def infer_sample_id_from_path(p: str) -> str:
    parent = os.path.basename(os.path.dirname(p))
    if not VALID_ID.match(parent):
        stem = os.path.splitext(os.path.basename(p))[0]
        parent = re.sub(r'_report$', '', stem)
    return parent

def normalise_id_to_sample_id(df: pd.DataFrame, want: str = 'Sample_ID') -> pd.DataFrame:
    if want in df.columns:
        return df
    for cand in ['sample_id','Sample_id','sample','Sample','isolate','Isolate','id','ID']:
        if cand in df.columns:
            return df.rename(columns={cand: want})
    return df

def read_qc_stack(qc_glob: str) -> pd.DataFrame:
    paths = sorted(glob.glob(qc_glob))
    if not paths:
        sys.exit(f"[combine] No QC files matched: {qc_glob}")
    dfs = []
    for p in paths:
        df = pd.read_csv(p, dtype=str)
        df = normalise_id_to_sample_id(df, 'Sample_ID')
        if 'Sample_ID' not in df.columns:
            df.insert(0, 'Sample_ID', infer_sample_id_from_path(p))
        df.replace("", pd.NA, inplace=True)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)

def read_typer(path_or_none: str) -> pd.DataFrame:
    if path_or_none == 'NONE':
        return pd.DataFrame(columns=['Sample_ID'])
    try:
        df = pd.read_csv(path_or_none, sep='\t', dtype=str)
    except Exception:
        df = pd.read_csv(path_or_none, sep=',', dtype=str)
    df = normalise_id_to_sample_id(df, 'Sample_ID')
    if 'Sample_ID' not in df.columns:
        sys.exit("[combine] Typer table missing a sample id column (e.g. Sample_ID)")
    df.replace("", pd.NA, inplace=True)
    return df

def ensure_columns(df: pd.DataFrame, cols: List[str]) -> None:
    for c in cols:
        if c not in df.columns:
            df[c] = pd.NA

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: generate_overall_report.py <qc_glob> <typer_path_or_NONE> <output_csv>")

    qc_glob, typer_path, out_csv = sys.argv[1], sys.argv[2], sys.argv[3]

    # 1) Load data
    qc_all   = read_qc_stack(qc_glob)
    typer_df = read_typer(typer_path)

    # 2) Build column order: fixed QC block + ALL typer columns (in typer file order, minus Sample_ID)
    typer_cols_in_order = [c for c in list(typer_df.columns) if c != 'Sample_ID']
    desired_order = QC_FIXED + typer_cols_in_order

    # 3) Outer-join on Sample_ID
    merged = pd.merge(qc_all, typer_df, on='Sample_ID', how='outer')

    # 4) Ensure fixed QC columns exist; select & order
    ensure_columns(merged, QC_FIXED)
    # Keep all columns: first QC_FIXED, then typer_cols_in_order, then any extras not yet included
    extras = [c for c in merged.columns if c not in set(['Sample_ID'] + QC_FIXED + typer_cols_in_order)]
    final_cols = ['Sample_ID'] + QC_FIXED[1:] + typer_cols_in_order + extras  # Sample_ID already first in QC_FIXED
    # Guarantee uniqueness and existence
    seen, ordered = set(), []
    for c in final_cols:
        if c in merged.columns and c not in seen:
            seen.add(c); ordered.append(c)

    out = merged.reindex(columns=ordered).sort_values('Sample_ID')

    # 5)
    #    - For PASS rows, missing in-silico fields => "MODULE FAILURE"
    if 'Overall_QC' in out.columns:
    mask = (out['Overall_QC'] == 'PASS')
    skip_cols = ['23S1_SNP', '23S3_SNP', 'gyrA_SNP', 'parC_SNP', 'typer_pipeline_version']
    fill_cols = [c for c in out.columns if c not in skip_cols]
    out.loc[mask, fill_cols] = out.loc[mask, fill_cols].fillna('MODULE FAILURE')

    #    - All remaining empties => "NA"
    out = out.fillna('NA')

    # 6) Write
    out.to_csv(out_csv, index=False)
    print(f"[combine] Wrote {out_csv} with {len(out)} rows × {len(out.columns)} cols")

if __name__ == "__main__":
    main()