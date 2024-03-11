#!/usr/bin/env python3

import pandas as pd


def main(fnames_coverage, samples, fname_all_coverage):

    tmp = []

    for file, sample in zip(fnames_coverage, samples):
        df = pd.read_csv(file, sep='\t')
        df['coverage'] = df[sample]

        tmp.appned(df[['pos', 'coverage']])

    pd.concat(tmp).to_csv(fname_all_coverage)


if __name__ == "__main__":
    main(
        snakemake.input.fnames_coverage,
        snakemake.params.samples,
        snakemake.output.fname_all_coverage,
    )
