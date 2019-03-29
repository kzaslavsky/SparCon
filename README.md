# SparCon

Code and data for analyses of sparse connectivity (SparCon) assays on stem cell derived neurons and other data (e.g., RNASeq) in Zaslavsky, Zhang et al. 2019.

SparCon assays are an extension of the co-culture approach first advanced in Shcheglovitov et al. (2013) intended to solve the problem of cell culture heterogeneity among stem cell-derived neurons that masks connectivity phenotypes. Specifically, SparCon assays use a consistent lawn of unlabaled neurons into which sparse numbers of differentially fluorescently labeled neurons are seeded. This enables tight control of the synaptogenic environment, as well as within-assay (e.g., within a single-well) normalization of assay (e.g., mutant) to benchmark (e.g., control) neurons, permitting comparisons across assays in separate wells, batches of wells, and perhaps different laboratories. Importantly, SparCon assays and within-well normalization significantly increase statistical power to detect differences in connectivity (e.g., synapse number, sEPSC frequency).

In addition to reproducing the plots shown in the main publication, this collection of files also contans generic code to allow users to analyze their own data by performing within-well normalization, Anderson-Darling k-samples tests, and generation of relevant plots.

We apologize for the messiness of the code. Any feedback and suggestions to improve the code are greatly appreciated!

=========

The provided R files are a collection of scripts and datasets used for analysis in the 
manuscript. Once run, they will re-create plots used in the manuscript.
Sample output is provided in folders.
Multiple tables of statistical analyses of the data are also provided in the output 
folders as .txt files

There are 5 folders:
'SparCon' contains files to perform analyses on SparCon data in the paper
'GENERIC_SPARCON_ANALYSIS' contains files to perform analyses on your own data.
'RNASeq' contains files to generate RNASeq analyses and figures in the paper.
'MiniScreen' contains files analyzing screen with chemical modulators in SparCon.
'Dendrite Extension' contains files analyzing the dendrite extension experiment

####### SparCon Folder ####### 

=== SCRIPTS ===

SparCon_analysis_v1.1.R 	- main analysis file. Will perform within-well normalization 
							for Synapse Number, Dendrite Length, sEPSC amplitude, 
							sEPSC frequency and produce the appropriate plots.
							- Requires: SparCon_fn.R

Distributions.R 			- fits distributions to datasets and produces corresponding 
							plots
							- Requires: SparCon_fn.R, pre-processing of data using 
							SparCon_analysis_v1.1.R

SparCon_instinsic_ephys.R 	- plots of electrophysiology data of intrinsic membrane 
							properties
							- Requires: SparCon_fn.R

Power_Sims.R 				- performs power simulations and outputs corresponding plots 
							and tables
							- Requires: SparCon_fn.R, Power_Sim_Fn.R, pre-processing of 
							data using SparCon_analysis_v1.1.R

Sholl_Analysis.R 			- parses Sholl dataset, performs within-well normalization, 
							outputs plots and performs statistical analysis
							- Requires: SparCon_fn.R

SparCon_fn.R 				- contains functions for within-well normalization, 
							statistical analyses, plotting

Power_Sim_Fn.R 				- contains functions for power simulations by sample size and
							fold change

=== DATASETS ===

SC_ICC_DB.txt 				- Imaging dataset containing Synapse Number and	Dendrite 
							Length data for 601 neurons. SparCon_analysis_v1.0.R, 
							Distributions.R, and Power_Sims.R use it. 

SC_ICC_DB_norm.txt 			- Within-well normalized dataset derived from
							SC_ICC_DB.txt

SC_ephys_db.txt				- Electrophysiology dataset containing sEPSC frequency and
							sEPSC amplitude data for 680 neurons. 
							SparCon_analysis_v1.0.R, Distributions.R, and Power_Sims.R 
							use it.
														
SC_ephys_db_norm.txt 		- Within-well normalized dataset derived from
							SC_ephys_db.txt
							
shollDB_unnorm.txt			- Raw dataset of Sholl analysis on 601 neurons. 
							Sholl_normalize.R parses it, performs within-well
							normalization, and outputs the Sholl analysis plots.

shollDB_norm.txt			- Within-well normalized dataset derived from
							shollDB_unnorm.txt

intrinsic_ephys_db.txt  	- electrophysiology dataset containing recordings of 
							intrinsic membrane properties. 
							Used by SparCon_instrinsic_ephys.R



####### RNASeq Folder ####### 

=== SCRIPTS ===

Deseq2_SHANK2.R          	- Differential expression analysis, prepares files for GSEA
							Generates output files used by other scripts

Activity_and_ASD_Modules.R  - Outputs log2 fold changes in activity-dependent transcripts
							and enrichment of ASD modules among DEGs

Volcano_plots.R      	    - generation of volcano plots and bargraphs of top and bottom 
							5 DEGs in each gene set for Fig 4

Heatmap_SHANK2.R            - generation of heatmap for Fig 4A

RNASeq_FN.R, SparCon_fn.R, GOfnKZ.R contain functions used by above files

=== DATASETS ===

simple_coding_counts_4wk.txt	- counts of aligned transcripts at 4wk

simple_coding_counts_9wk.txt    - counts of aligned transcripts at 9wk

Pruunslid_Bic_Stim_7wk.txt      - activity-dependent transcripts from Pruunslid et al.

Neuron_Markers.txt              - manually-curated markers of cortical layers and synapses

FMRP.csv                        - FMRP targets from Darnell et al.

ASD_Modules.txt                 - Modules from WGCNA analyses by Parikshak et al. 
								and Willsey et al.

geneset_list.csv                - genesets and GOIDs to facilitate plotting Fig 4


####### MiniScreen Folder ####### 

=== SCRIPTS ===

Miniscreen_plots_script.R		- script to generate figures for the MiniScreen

=== DATASETS ===

Compiled_miniscreen_data_raw.csv	- raw miniscreen data for 1246 neurons

Mutant_normalized_full.csv and Mutant_normalized_stats.csv - data and statistics 
normalized to R841X DMSO

WT_normalized_data.csv and WT_normalized_stats.csv - data and statistics normalized to 
R841X-C DMSO

Within_genotype_normalized_data.csv and Within_genotype_normalized_stats.csv - data
and statistics normalized within each genotype

####### Dendrite Extension Folder ####### 

=== SCRIPTS ===

dend_extension_exp.R         	- script to generate figures from the dendrite extension 
								experiment.

=== DATASETS ===

Dend_extend_baseline.txt		- dataset containing dendrite length at baseline

Dend_extend_normalize.txt       - dataset containing dendrite lengths at 500 min 
                                normalized to dendrite length at baseline for each neuron

####### GENERIC_SPARCON_ANALYSIS Folder ####### 

=== SCRIPTS ===

Generic_SparCon.R            	- Change variable names & plotting parameters as desired. 
								Enter your own data as in the example text files.

Generic_SparConFN.R				- functions to normalize within-well, output stats and 
								plots.
