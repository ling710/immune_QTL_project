# immune_QTL_project

## Description of the project
We applied two-sample Mendelian randomization (MR) and colocalization to identify candidate genes in different immune cell types with putative causal roles on cancer outcomes. We expect top findings could have associations with immune checkpoints or genes in the immune response pathway or immune modulators and illustrate it in network plot. We expect to classify people into different groups based on our top findings, to detect characteristics and survival status and specific genes’ expression levels in different groups.

## Description of the files in the repository
We have five folders: scripts for code used in the project
                      raw_data for original data links
                      processed_data for any processed data for MR analysis, network analysis and plots
                      results for key tables and plots in the project
                      docs for any issue discussed

## how to run the scripts in the project and get the corresponding results 
1.clone the repository, download the raw data according to the links in the raw_data folder.
2.run DICE_data_preprocessing.R for data processing
3.run DICE_MR.R for 2SMR analysis, steiger filtering,colocalization and ldcheck
4.run DICE_plot.R for MR forest plots
5.run DICE_pwcoco.sh for PWCOCO analysis
6.run network.R for network analysis and plots
oneK1K workflow is similar to DICE






