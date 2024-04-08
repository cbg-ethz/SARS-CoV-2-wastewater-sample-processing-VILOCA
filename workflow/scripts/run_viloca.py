import subprocess
from pathlib import Path
from os import listdir
from os.path import isfile, join
import pandas as pd


def main(fname_bam, fname_reference, fname_insert_bed, fname_results_snv, dname_work, fname_bad_samples, sample, threads):

    bad_samples = pd.read_csv(fname_bad_samples)['sample'].values.tolist()

    if sample in bad_samples:
        with open(fname_results_snv.resolve(), 'w') as f:
            f.write('#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO')
        #f.close()

    else:

        alpha = 0.000001
        n_max_haplotypes = 100
        n_mfa_starts = 1

        dname_work.mkdir(parents=True, exist_ok=True)

        subprocess.run(
            [
                "viloca",
                "run",
                "-b",
                fname_bam.resolve(),
                "-f",
                fname_reference.resolve(),
                "--mode",
                "use_quality_scores",
                "--alpha",
                str(alpha),
                "--n_max_haplotypes",
                str(n_max_haplotypes),
                "--n_mfa_starts",
                str(n_mfa_starts),
                #'--exclude_non_var_pos_threshold',
                #str(0.001),
                "-z",
                fname_insert_bed.resolve(), # amplicon mode doesn't work because we couldn't merge reverse and forward reads
                #"-w",
                #str(249),
                #"--region",
                #"NC_045512.2:10099-23327",
                "--min_windows_coverage",
                "1",
                "--threads",
                str(threads)
            ],
            cwd=dname_work,
        )

        (dname_work / "snv" / "SNVs_0.010000_final.vcf").rename(fname_results_snv)


if __name__ == "__main__":
    main(
        Path(snakemake.input.fname_bam),
        Path(snakemake.input.fname_reference),
        Path(snakemake.input.fname_insert_bed),
        Path(snakemake.output.fname_vcf),
        Path(snakemake.output.dname_work),
        snakemake.input.fname_bad_samples,
        snakemake.params.sample,
        snakemake.threads
    )
