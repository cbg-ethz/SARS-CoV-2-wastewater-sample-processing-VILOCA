"""
Input:
    -- snvs.vcf outputted by VILOCA
    -- yaml with mutation definitions

Output:
    -- per sample csv file with the mutation calls and their respective frequency
"""
import yaml
import pandas as pd
from fuc import pyvcf

def parse_yaml(fname_yaml):

    yaml_file = open(fname_yaml)
    parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)
    dict_mut =parsed_yaml_file.get('mut')

    return dict_mut





def main(fname_vcf, fname_csv, fname_yaml):

    # read mutations
    dict_mut = parse_yaml(fname_yaml)
    positions_of_interest = [int(x) for x in list(dict_mut.keys())]

    if fname_vcf.split(".")[-1]=="vcf":
        # vcf into dataframe
        df_vcf = pyvcf.VcfFrame.from_file(fname_vcf).df
        df_mut = df_vcf[df_vcf['POS'].isin(positions_of_interest)]

    elif  fname_vcf.split(".")[-1]=="tsv":
        df_vcf = pd.read_csv(fname_vcf, sep='\t')
        print(df_vcf.columns())
        df_mut = df_vcf[df_vcf['Pos'].isin(positions_of_interest)]

    df_mut.to_csv(fname_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fname_snvs_vcf,
        snakemake.output.fname_mutations,
        snakemake.params.fname_mutation_definitions,
    )
