"""
Input:
    -- csv files of all mutations
    -- yaml with mutation definitions

Output:
    -- per sample csv file with the mutation calls and their respective frequency
"""
import yaml
import pandas as pd

def parse_yaml(fname_yaml):

    yaml_file = open(fname_yaml)
    parsed_yaml_file = yaml.load(yaml_file, Loader=yaml.FullLoader)
    dict_mut =parsed_yaml_file.get('mut')

    return dict_mut


def main(fname_muts, fname_csv, fname_yaml, all_samples):

    # read mutations
    dict_mut = parse_yaml(fname_yaml)
    positions_of_interest = [int(x) for x in list(dict_mut.keys())]

    df_muts = pd.read_csv(fname_muts)

    # filter for positions of interest
    df_muts = df_muts[df_muts['POS'].isin(positions_of_interest)]

    # for samples that don't have positins of interest we need to add an empty lind
    for sample_name in all_samples:
        if sample_name not in df_muts['sample'].unique():
            # create empty row
            df_muts = df_muts.append({'sample':sample_name}, ignore_index=True)
            #df_muts = pd.concat([df_muts, pd.DataFrame({"sample": sample_name})])

    df_muts.to_csv(fname_csv)

if __name__ == "__main__":
    main(
        snakemake.input.fname_all_mutations,
        snakemake.output.fname_mutations,
        snakemake.params.fname_mutation_definitions,
        snakemake.params.all_samples,
    )
