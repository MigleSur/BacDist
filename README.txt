----- General information -----

All the code is in the Snakefile and is written in snakemake.

The pipeline takes a number of WGS bacterial genome fastq files and outputs:
	(1) the VCF file with unique(minimum one sample doesn't have the variant), high quality variants
	(2) SNP-distance matrix based on the filtered VCF file
	(3) phylogenetic tree (if more than 3 samples were used for the analysis) based on the filtered VCF file and created using RAxML

----- Setting up the config file -----

In order for the pipeline to run, a configuration file is needed. A configuration file requires 4 fields to be present: 
	* sample_dir - Directory in which all fastq files which need to be analyzed are present.
	* output_dir - Output directory where all the output files will be created.
	* ref - Reference genome in gbk format. 
	* name - Chosen name for the analysis. The name is used as a prefix in output and intermediate files.

Here is the example configuration file:
sample_dir: "/home/projects/cu_10047/people/misur/PhD_projects/Achromobacter/fastq_files/AX01_input"
output_dir: "/home/projects/cu_10047/people/misur/PhD_projects/Achromobacter/variant_analysis/AX01_output"
ref: "/home/projects/cu_10047/people/misur/PhD_projects/Achromobacter/reference_genomes/GCF_001457475.1_NCTC10807_genomic.gbk"
name: "AX01"

----- Running the pipeline -----

In order to run the pipeline anaconda3/4.0.0 module has to be loaded. Snakemake is started from its directory (/home/projects/cu_10047/people/misur/PhD_projects/scripts/multiple_bacteria_SNP_analysis/) 
	
	-j option allows to choose the number of threads (1-28) used for the analysis (default:1)
	--configfile option allows to chose the configuration file for the analysis
	
Here is the example code for running snakemake:

	module load anaconda3/4.0.0
	snakemake -j 10 --configfile config.yaml 

----- Rerunning the pipeline -----

Pipeline can be rerun when new samples are added, old samples are deleted, etc. 

First, all the new samples have to be added to the sample directory.
Second, the following directories have to be deleted:
	* temp
	* output_files
	*vcf_calls 
After completing the steps a pipeline can be rerun with different samples.
