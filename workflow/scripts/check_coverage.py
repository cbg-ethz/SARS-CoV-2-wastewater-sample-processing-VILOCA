#!/usr/bin/env python3

import pandas as pd


def main(fnames_coverage, region_of_interest, samples, fname_good_samples, fname_bad_cov_samples):

    bad = []
    good = []

    start_region = int(region_of_interest[0])
    end_region = int(region_of_interest[1])
    region_list = list(range(start_region, end_region))

    for file, sample in zip(fnames_coverage, samples):
        df = pd.read_csv(file, , sep='\t')
        df = df.loc[df['pos'].isin(region_list)]

        if df[sample].mean() < 100:
            bad.append(sample)
        else:
            # pass coverage test
            good.append(sample)

    df_bad = pd.DataFrame()
    df_bad['sample']=bad
    df_bad.to_csv(fname_bad_cov_samples)

    df_good = pd.DataFrame()
    df_good['sample']=good
    df_good.to_csv(fname_good_samples)


if __name__ == "__main__":
    main(
        snakemake.input.fnames_coverage,
        snakemake.params.region_of_interest,
        snakemake.params.samples,
        snakemake.output.fname_samples,
        snakemake.output.fname_bad_cov_samples,
    )
