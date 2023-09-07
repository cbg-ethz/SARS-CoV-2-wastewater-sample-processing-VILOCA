#!/usr/bin/env python3

def main(input_files, output_file):

    with open(output_file, 'w') as outfile:
        for fname in input_files:
            with open(fname) as infile:
                outfile.write(infile.read())



if __name__ == "__main__":
    main(
        snakemake.input.haplos,
        snakemake.output.fname,
    )
