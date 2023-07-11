# Wastewater sample processing

This workflow processes wastewater samples for the detection of drug resistance mutations.

Input:
- aligned samples (alignment with V-pipe)
- positions of interest: positions where drug resistance mutations are expected
- list of all samples to be processed

Output:
- csv file of all mutations at positions of interest
--> This file can be used in the visualisation notebook to create a heatmap

To run this workflow adapt file paths in config and sample names. Run the pipeline on Euler with run_workflow.sh.

Note, at the moment the workflow is SARS-CoV-2 specific, however, reference files and gene annotations can be easily swapped. 
