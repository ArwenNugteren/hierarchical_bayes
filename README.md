This repository contains the code and data for the paper 'Hierarchical Bayes meets hierarchical forecasting: A flexible framework for level-focused forecasts'. 

Included are the files to: 
- generate the synthetic data (both coherent and incoherent versions);
- run the simulation study models;
- run the tourism data models;
- compute accuracy metrics for individual trials, and those to collate the results into the tables and plots used in the paper.


## Reproducing the paper
The file `required_packages.R` installs the packages required for all files in this repository. Packages and version used in writing this paper: 
```
> sessionInfo()
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] MASS_7.3-65         scoringRules_1.1.3  tidyr_1.3.1        
 [4] fable_0.4.1         fabletools_0.5.0    tsutils_0.9.4      
 [7] cowplot_1.2.0       patchwork_1.3.2     ggplot2_3.5.2      
[10] rstan_2.32.7        StanHeaders_2.32.10 stringr_1.5.1      
[13] dplyr_1.1.4
```

### Simulation study
All files located in the `simulation_study` folder. The folder `data` contains the generated data and true parameters. The file `synthetic_hierarchy_hierarchical_params.R` generates the raw data, and the file `incoherent_synthetic_hierarchy.R` generates the softly coherent data. The folder `prior_predictive_checks` contains the code and images in the appendix for the misspecified and wellspecified independent and hierarchical Bayesian models. 

Models are run in the files: 
- `hierarchy_independent_misspecified.R` and `hierarchy_misspecified.R` run the independent and hierarchical Bayesian trials on the misspecified data with MinT reconciliation. `hierarchy_misspecified_coherency.R` runs the coherency and level-weighted trials by editing the values of `sigma_coh`, `lambda_top` and `lambda_bot`.
- `hierarchy_wellspecified_joint_independent.R` and `hierarchy_wellspecified_joint.R` run the independent and hierarchical Bayesian trials on the misspecified data with MinT reconciliation. `hierarchy_wellspecified_joint_coherency.R` runs the coherency and level-weighted trials by editing the values of `sigma_coh`, `lambda_top` and `lambda_bot`.
For the trials with varying window size and incoherent data, edit the set-ups in Section 2 of the files to uncomment the dataset you wish to use and to slice it and edit the training window size. While individual windows can be run on a laptop computer in under 2 minutes, a HPC system was used to run 50 batches of 10 windows each to save time over the full 500 window rolling evaluation. A sample PBS job script is shown. Each batch on the HPC ran around 5 minutes, or up to 15 minutes for coherency penalised trials. 5GB of memory per batch was found to be sufficient. 10 CPUS were used for the 10 chains run. 

Accuracy results: 
- Tables in the paper are obtained by running `tables_misspecified.R` and `tables_wellspecified.R` which save to the folders `tables_misspecified` and `tables_wellspecified` respectively, under the file names `overall.table.csv` (upper portion of Tables 3 and A.9), `results_level_table.csv` (lower portion of Tables 3 and A.9) and `MAI.csv` (Tables 2 and A.10). Tables 1 and A.8 are obtained from the the same files (reading in different accuracy folders corresponding to the appropriate trials) and percentage improvements were calculated in Excel from the output `overall_table.csv`. 
- `comparison_misspecified_incoherent.R` and `comparison_wellspecified_incoherent.R` summarise the results for the incoherent data and save data to the folders `results_misspecified_incoherent` and `results_wellspecified_incoherent`. The paper tables 4, 5, A.11 and A.12 are obtained by combining these with the results on the coherent data calculated previously.
The file `accuracy_synthetic.R` can be used to obtain the accuracy results for a single model as opposed to comparing multiple. 




### Australian domestic tourism example
Files are in the `australian_tourism` folder. This dataset is originally obtained from the data in the paper of Wickramasuriya et al. (2019). 

Data, raw and processed, is found in the data folder. The file `data_prep.R` cleans up the data and cleans it into an appropriate format, saved as `tourism.Rdata`, which is called by all other scripts. Prior predictive checks are in the `prior_predictive_checks` folder. The scripts `prior_predictive_check_tourism_independent.R` and `prior_predictive_check_tourism_hierarchical.R` run the two prior predictive checks and saves the output images to the same folder. Both the data cleaning and prior predictive checks can be run locally on a regular laptop computer. 

Models are run in the files below. Each file saves fitted models and accuracy metrics locally (these are too large to upload to github).
- `purposes_independent_MinT_rolling`: runs the Independent Bayes models and the MinT reconciled model.
- `purposes_hierarchy_MinT_rolling`: runs the Hierarchical Bayes model and MinT reconciled version. 
- `purposes_hierarchy_coherency_rolling`: runs coherency weighted trials. The value of `sigma_coh` can be edited to give the softer or tighter coherency penalties.
- `purposes_hierarchy_coherency_rolling`: runs coherency weighted trials. The values of `lambda_total` can be edited to `n_rp` for the top-level weighting and `lambda_pp` can be edited to `n_rp/n_pp` for the purpose level weighting.

An example PBS script, `tourism_job.pbs` for running batch jobs on HPC is provided. This script creates folders with the job name. Packages must be installed in your HPC environment first (instructions are not provided as HPC systems differ). The independent and hierarchical trials ran within 4 hours per batch on HPC, while the coherency penalised trials took up to 6 hours due to the more complex geometry of the posterior induced by penalisation. 96GB of memory was found to be sufficient with some leeway. 10 CPUS were used for the 10 chains of the MCMC sampler. 

Results are obtained by running the file `paper_tables.R`. This creates and saves results to the folder `paper_tables`. The file `overall_results.csv` contains the RMSSE, ES and 95% coverage at the top of Table 6, while the file `level_results.csv` contains the RMSE, CRPS and 95% coverage by level in the lower part of Table 6. The file `MAI.csv` contains the results in Table 7. The tables in the appendix are calculated separately in the file `appendix_tables.R`, save to the folder `appendix_tables`, and contains results by horizon in the files `RMSSE.csv`, `ES.csv`, `RMSE.csv` and `CRPS.csv`. The file `paper_plots.R` creates the plots in the tourism data example, and saves them to a folder called `paper_plots`. This script creates a pdf per metric, containing the Nemenyi results for ES and RMSSE (Figure 3), followed by the relative line plots by level for RMSE and CRPS (Figures 4 and 5). In all of these files, the user will need to point towards the appropriate file paths for the accuracy folders (either copied over locally, or in a HPC drive). 

The file `accuracy_purposes_rolling.R` runs accuracy metrics for a single trial, and the file `post_hoc_es` contains some additional files that were used due to an earlier error in energy score calculation in the model files, so energy score results were re-run post-hoc from the saved models. This error has now been rectified, but the files are retained for completeness. 


Corresponding author for queries on this repository is Arwen Nugteren, email: a.nugteren@qut.edu.au
This repository last edited 2nd October 2026. 

### References
Wickramasuriya, S.L., Athanasopoulos, G., Hyndman, R.J., 2019. Optimal forecast reconciliation for hierarchical and grouped time series through trace minimization. Journal of the American Statistical Association 114, 804–819. doi:10.1080/01621459.2018.1448825

