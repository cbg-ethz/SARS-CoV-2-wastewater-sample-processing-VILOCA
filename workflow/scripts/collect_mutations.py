#!/usr/bin/env python3

"""
Script collecting the mutations of interest from all samples.
"""

import pandas as pd

def concat_files(in_files_list, out_fname, params):
    tmp = []

    for file, param in zip(in_files_list, params):
        df = pd.read_csv(file)
        df["sample"] = param
        tmp.append(df)

    pd.concat(tmp).to_csv(out_fname)


def main(csv_list, fname_result_csv, params):

    concat_files(csv_list, fname_result_csv, params)


if __name__ == "__main__":
    main(
        snakemake.input.csv_list,
        snakemake.output.fname_result_csv,
        params.params,
        )
