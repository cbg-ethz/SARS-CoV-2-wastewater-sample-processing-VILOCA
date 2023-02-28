import pandas as pd
import subprocess
import os


def subsample(fname_bam_in, fname_bam_out , subsample_proportion):
    """
    from samtools view:
    --subsample: Output only a proportion of the input alignments, as specified
    by 0.0 ≤ FLOAT ≤ 1.0, which gives the fraction of templates/pairs to be kept.
    This subsampling acts in the same way on all of the alignment records in the
    same template or read pair, so it never keeps a read but not its mate.
    """

    subprocess.run(
        [
            "samtools",
            "view",
            "-s",
            str(subsample_proportion),
            "-b",
            str(fname_bam_in),
            "-o",
            str(fname_bam_out),
        ],
    )

def get_subsample_proportion(fname_coverage, sample, coverage_threshold=100000):
    """
    viloca threshold is coverage_threshold=1000000
    """

    df = pd.read_csv(fname_coverage, sep='\t')
    max_coverage = df[sample].max()

    if max_coverage <= coverage_threshold:
        return 1
    else:
        return coverage_threshold/max_coverage


def main(fname_bam_in, fname_bam_out, fname_coverage, sample, coverage_threshold):

    subsample_proportion = get_subsample_proportion(fname_coverage, sample, coverage_threshold)
    if subsample_proportion==1:
        os.rename(fname_bam_in, fname_bam_out)
    else:
        subsample(fname_bam_in,fname_bam_out,subsample_proportion)


if __name__ == "__main__":
    main(
        snakemake.input.fname_bam,
        snakemake.output.fname_bam,
        snakemake.input.fname_coverage,
        snakemake.params.sample,
        snakemake.params.coverage_threshold,
    )
