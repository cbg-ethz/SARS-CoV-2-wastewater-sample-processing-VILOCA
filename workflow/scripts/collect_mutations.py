#!/usr/bin/env python3

"""
Script collecting the mutations of interest from all samples.
"""

import pandas as pd

def concat_files(in_files_list, params):
    tmp = []

    for file, param in zip(in_files_list, params):
        df = pd.read_csv(file)
        df["sample"] = param
        tmp.append(df)

    return pd.concat(tmp)

def split_INFO(df_temp):
    # get parameters
    df_temp['INFO_list'] = df_temp['INFO'].str.split(";")
    parameter_list = [ parameter.split('=')[0] for parameter in df_temp['INFO_list'][0]]
    for param in parameter_list:
        df_temp[param]=' '

    for iter_row, row in df_temp.iterrows():
        for param in parameter_list:
            df_temp.loc[iter_row,param]=row['INFO'].split(param+"=")[1].split(";")[0]

    return df_temp


def main(csv_list, fname_result_csv, params):

    df = concat_files(csv_list, params)
    df = split_INFO(df)
    df.to_csv(fname_result_csv)


if __name__ == "__main__":
    main(
        snakemake.input.csv_list,
        snakemake.output.fname_result_csv,
        snakemake.params.params,
        )
