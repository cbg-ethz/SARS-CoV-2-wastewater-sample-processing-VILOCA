import pysam
from pathlib import Path
import pandas as pd

def count_read(bam_file):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Initialize a count variable
    read_count = 0

    # Iterate over each read in the BAM file and increment the count
    for read in bam:
        read_count += 1

    # Close the BAM file
    bam.close()
    return read_count


def main(fnames_bam, fnames_bam_fused, fnames_bam_nonfused, fname_csv)
    collect_info = {}

    for fname_bam, fname_bam_fused, fname_bam_nonfused in zip(fnames_bam, fnames_bam_fused, fnames_bam_nonfused):

        collect_info.update({
            'sample': str(fname_bam).split("results/")[1].split("/alignment")[0],
            'read_count_before': count_read(fname_bam),
            'read_count_fused': count_read(fname_bam_fused),
            'read_count_nonfused': count_read(fname_bam_nonfused),
        })

    pd.DataFrame(collect_info).to_csv(fname_csv)




if __name__ == "__main__":
    main(
        Path(snakemake.input.fnames_bam),
        Path(snakemake.input.fnames_bam_fused),
        Path(snakemake.input.fnames_bam_nonfused),
        Path(snakemake.output.fname_csv),
    )
