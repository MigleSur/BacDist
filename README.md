# BacDist
Snakemake pipeline for bacterial SNP distance and phylogeny analysis

## General information

All the code is in the Snakefile and is written in snakemake.

The pipeline takes a number of WGS bacterial genome fastq files and outputs:

1. VCF file with unique (minimum one sample doesn't have the variant), high quality variants
	
2. SNP-distance matrix based on the filtered VCF file
	
3. Phylogenetic tree (if more than 3 samples were used for the analysis) based on the filtered VCF file and created using RAxML

## Required software

Before running the pipeline, make sure that the following programs are installed and added to the path:

[freebayes>=1.1.0]()

[vcflib>=1.0.0-rc2]()

[vcftools>=0.1.16]()

[snpeff>=4.3r]()

[prokka>=1.12]()

[minimap2>=2.6]()

[seqtk>=1.0]()

[snp-sites>=2.4.0]()

[emboss>=6.6.0]()

[bcftools>=1.9]()

[snippy>=4.1.0](https://github.com/tseemann/snippy)

[vt>=0.5772]()

[bcftools>=1.9]()

[samtools>=1.9]()

[vcftools>=0.1.16]()

[raxml>=8.2.11]()


### Disclaimer

[Snippy4](https://github.com/tseemann/snippy) doesn't work with python3. Python3 should be disabled at that step and Python2 should be available.

## Setting up the config.yaml file 

In order for the pipeline to run, a configuration file is needed. A configuration file requires 4 fields to be present: 
* sample_dir - Directory in which all fastq files which need to be analyzed are present.
* output_dir - Output directory where all the output files will be created.
* ref - Reference genome in gbk format. 
* name - Chosen name for the analysis. The name is used as a prefix in output and intermediate files.

Here is the example configuration file:
```
sample_dir: "/home/project_name/fastq"
output_dir: "/home/project_name/AX01"
ref: "/home/reference_genomes/GCF_001457475.1_NCTC10807_genomic.gbk"
name: "AX01"
```

## Running the pipeline 

In order to run the pipeline anaconda3 (version 4.0.0) has to be available. Snakemake is started from its directory directory:
    
```
-j option allows to choose the number of threads (1-28) used for the analysis (default:1)
--configfile option allows to chose the configuration file for the analysis
```

Here is the example code for running snakemake:

```bash
	snakemake -j 10 --configfile config.yaml 
```

## Rerunning the pipeline 

Pipeline can be rerun when new samples are added, old samples are deleted, etc. 

First, all the new samples have to be added to the sample directory.

Second, the following directories have to be deleted:
* temp
* output_files
* vcf_calls 

After completing the steps a pipeline can be rerun with different samples.

## Author

Migle Gabrielaite | migle.gabrielaite@gmail.com
