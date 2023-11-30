#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import glob
import subprocess

def unzip(fname):
    subprocess.run(
            [
                "gunzip",
                str(fname),
            ],
        )


def main(dname_haplotypes, sample_name, fname_out_haplos):

    for fname in glob.glob(f"{dname_haplotypes}*.reads-support.fas.gz"):
        unzip(fname)

    all_haplotype_files = [file for file in glob.glob(f"{dname_haplotypes}*.reads-support.fas")]

    # define output file
    with open(fname_out_haplos, 'w') as out_f:
        my_seqs = []

        # iterate over all haplotype files
<<<<<<< HEAD
        for fname_haplotype_region in all_haplotype_files:
=======
        for fname_haplotype_region in fnames_haplotypes:
>>>>>>> 2fe8abf0436151a3940e8b84cfdfff39bf636c8e
            region = fname_haplotype_region.split(".reads-support.fas")[0].split("w-")[1]
            for record in SeqIO.parse(fname_haplotype_region, 'fasta'):
                    new_id = sample_name+'-'+region + '-' +str(record.id)
                    my_seqs.append(SeqRecord(Seq(record.seq), id = new_id))

        # write all sequences to file
        SeqIO.write(my_seqs, out_f, "fasta")



if __name__ == "__main__":
    main(
        snakemake.input.dname_haplotypes,
        snakemake.params.sample,
        snakemake.output.fname_haplotypes,
    )
