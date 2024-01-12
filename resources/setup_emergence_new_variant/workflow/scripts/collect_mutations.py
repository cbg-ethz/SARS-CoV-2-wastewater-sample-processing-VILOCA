#!/usr/bin/env python3
"""
Script aggregating co-occurring mutations from all samples.
"""
import pandas as pd

def main(fnames_csv, fname_cooccurring_mutations_csv):

    tmp = []

    for fname_csv in fnames_csv:

        sample = fname_csv.split("results/")[1].split("/variant_calling/snv/")[0]

        df_tmp = pd.read_csv(fname_csv)
        df_tmp['sample'] = sample

        tmp.append(df_tmp)

    merged_csv = pd.concat( tmp )
    merged_csv.to_csv(fname_cooccurring_mutations_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fnames_csv,
        snakemake.output.fname_csv,
    )
