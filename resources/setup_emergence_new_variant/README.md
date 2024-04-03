This workflow produces Figure 6 in the Manuscript of VILOCA.

To reproduce the results and figure:  
1.) Clone the repository.
2.) Move into this directory: `cd resources/setup_emergence_new_variant`
3.) Install conda enviroments needed for the the workflow:`snakemake --conda-create-envs-only --use-conda -c1 --rerun-incomplete`
4.) Execute workflow. On slurm cluster the script `run_workflow.sh` can be used.
5.) When you received the two needed output files `results/all_cooccurring_mutations.csv` and `results/samples.all_coverage.csv`, figure can be generated with the notebook `workflow/notebooks/rise_of_ba.5.ipynb`
