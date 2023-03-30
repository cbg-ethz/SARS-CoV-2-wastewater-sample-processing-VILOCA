#!/usr/bin/env python3
"""
Collect annotated mutations from all samples.
"""
import pandas as pd
from fuc import pyvcf


def concat_files(in_files_list, params):
    tmp = []

    for file, param in zip(in_files_list, params):

        df = pyvcf.VcfFrame.from_file(file).df
        df["sample"] = param
        df = split_INFO(df)

        tmp.append(df)

    return pd.concat(tmp)

def split_INFO(df_temp):
    # get parameters
    df_temp['INFO_list'] = df_temp['INFO'].str.split(";")
    first_row_INFO_list = df_temp['INFO_list'].values.tolist()[0]
    parameter_list = [param_info.split('=')[0] for param_info in first_row_INFO_list]
    for param in parameter_list:
        df_temp[param]=' '

    for iter_row, row in df_temp.iterrows():
        for param in parameter_list:
            df_temp.loc[iter_row,param]=row['INFO'].split(param+"=")[1].split(";")[0]

    return df_temp


def main(vcf_list, fname_result_csv, params):

    df = concat_files(vcf_list, params)
    df= df.reset_index()

    df.to_csv(fname_result_csv)


if __name__ == "__main__":
    main(
        snakemake.input.vcf_list,
        snakemake.output.fname_result_csv,
        snakemake.params.params,
        )
