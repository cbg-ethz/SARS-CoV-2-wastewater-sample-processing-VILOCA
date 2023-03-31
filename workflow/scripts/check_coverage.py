#!/usr/bin/env python3

import pandas as pd


def main(fnames_coverage, samples, fname_good_samples, fname_bad_cov_samples):

    bad = []
    good = []

    for file, sample in zip(fnames_coverage, samples):
        df = pd.read_csv(file, sep='\t')

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
        snakemake.params.samples,
        snakemake.output.fname_samples,
        snakemake.output.fname_bad_cov_samples,
    )
