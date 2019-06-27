import sys
import argparse
import glob, os


if os.path.exists('config.yaml'):
        configfile: 'config.yaml'

input_dir=config['sample_dir']
reference=config['ref']
output=config['output_dir']


def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)



all_fastqs = glob.glob(input_dir+"/*R1*")
temp_fastqs = [fastq.rsplit('_R1',1)[0] for fastq in all_fastqs]
IDS = [fastq.replace(input_dir+'/', '') for fastq in temp_fastqs]


sample_no=len(IDS)
sample_no2=len(IDS)*2


name_with_ref=''.join((config["name"],"_incl_reference"))

# run all the jobs to produce the final files
rule all:
	input:
		final_vcf=expand("{outdir}/output_files/{name}_final.vcf", name=config["name"], outdir=config["output_dir"]),
		snp_dist=expand("{outdir}/output_files/{name}_snp_dist.tsv", outdir=config["output_dir"], name=config["name"]),
		raxml=expand("{outdir}/output_files/RAxML_bestTree.{name}", outdir=config["output_dir"], name=config["name"]),
		raxml_with_ref=expand("{outdir}/output_files/RAxML_bestTree.{name}", outdir=config["output_dir"], name=name_with_ref)
	message:
		"Pipeline complete"


def check_fastq_end():
	fastqs = glob.glob(input_dir+"/*R1*")
	fastqs += glob.glob(input_dir+"/*R2*")
	end = [fastq.rsplit('_R', 1)[1][1:] for fastq in fastqs]
	if len(set(end)) == 1:
		return(list(set(end))[0])
	else:
		raise Exception("ERROR file endings are different. Make sure that all the files have the same file endings.")
file_end=check_fastq_end()



# running snippy4 on all input samples
rule run_snippy:
	input:
		R1=expand("{input_dir}/{{sample}}_R1{end}", input_dir=config["sample_dir"], end=file_end),
		R2=expand("{input_dir}/{{sample}}_R2{end}", input_dir=config["sample_dir"], end=file_end)
	output:
		vcf=expand("{outdir}/raw_vcf_calls/{{sample}}/{{sample}}.vcf.gz", outdir=config["output_dir"]),
		bam=expand("{outdir}/raw_vcf_calls/{{sample}}/{{sample}}.bam", outdir=config["output_dir"]),
		ref=expand("{outdir}/raw_vcf_calls/{{sample}}/reference/ref.fa", outdir=config["output_dir"])
	params:	
		outdir=expand("{outdir}",outdir=config["output_dir"]),
		ref=config['ref'],
		prefix="{sample}",
		cpus=10,
		mapqual=50,
		mincov=10,
		minfrac=0.5,
		minqual=50
	log:
		expand("{outdir}/logs/snippy_{{sample}}.log",outdir=config["output_dir"])
	message:
		"Running snippy on all provided samples"
	shell:
		"""
                module unload anaconda3/4.0.0
		module load anaconda2/4.0.0
                module load tbl2asn/20190117
		module load parallel/20190122
		module load bwa/0.7.15
                module load samtools/1.9
		module load java/1.8.0
		module load jre/1.8.0		
		module load perl/5.24.0
                module load freebayes/1.1.0-50-g61527c5
                module load vcflib/1.0.0-rc2
		module load vcftools/0.1.16
                module load snpeff/4.3r
                module load prokka/1.12
                module load minimap2/2.6
                module load seqtk/1.0-r82-dirty
                module load snp-sites/2.4.0
                module load emboss/6.6.0
                module load bcftools/1.9
                module load snippy/4.1.0
		module load vt/0.5772
	
		(snippy --outdir {params.outdir}/raw_vcf_calls/{params.prefix} --ref {params.ref} --R1 {input.R1} --R2 {input.R2} --prefix {params.prefix} --cpus {params.cpus} --mapqual {params.mapqual} --mincov {params.mincov} --minfrac {params.minfrac} --minqual {params.minqual}  --force) 2> {log}
		"""


input_merge = [''.join((output,'/raw_vcf_calls/',file,'/',file,'.vcf.gz')) for file in IDS]



# Merge all variant files
rule merge_vcfs:
	input:
		calls=input_merge	
	output:
		expand("{outdir}/vcf_calls/merged_{name}_raw.vcf", name=config["name"], outdir=config["output_dir"])
	log:
		expand("{outdir}/logs/vcf_merge.log",outdir=config["output_dir"])
	message:
		"Merging all VCF files"
	shell:
		"""
		module load bcftools/1.9
		
		(bcftools merge -0 -o {output} {input.calls}) 2> {log}
		"""


# Filtering the variants in the merged file
rule filter_out_same_variants:
	input:
		expand("{outdir}/vcf_calls/merged_{name}_raw.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		expand("{outdir}/vcf_calls/{name}_diff.vcf", name=config["name"], outdir=config["output_dir"])
	params:
		filter="\"( GEN[*].GT='0/0' )\""
	message:
		"Filtering SNPs for all provided samples: Excluding variants shared by all isolates."
	shell:
		"""
		snpsift=/services/tools/snpeff/4.3r/snpEff/SnpSift.jar
		java -jar $snpsift filter {params.filter} {input} > {output}
	
		"""


# creating a bed file with all positions in the merged VCF file
rule create_all_position_bed:
	input:
		vcf=expand("{outdir}/vcf_calls/{name}_diff.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		bed=temp(expand("{outdir}/temp/{name}_temp.bed", name=config["name"], outdir=config["output_dir"]))
	shell:
		"""
		cat {input.vcf} | grep -v "#" | awk '{{print $1,$2-1,$2}}' > {output.bed}

		"""

# Extracting the low quality positions
rule extract_low_coverage_positions:
	input:
		bam=expand("{outdir}/raw_vcf_calls/{{sample}}/{{sample}}.bam", outdir=config["output_dir"]),
		bed=expand("{outdir}/temp/{name}_temp.bed", name=config["name"], outdir=config["output_dir"])
	output:
		low_qual=temp(expand("{outdir}/temp/{{sample}}_low_qual",  outdir=config["output_dir"]))
	shell:
		"""
		module load samtools/1.8
		samtools depth -b {input.bed} {input.bam} | awk '$3<10 {{print $1,$2-1,$2}}' | tr " " "\t" > {output.low_qual}
		"""

# Temp rule to avoid ambiguity
rule copy_input_to_output:
	input:
		vcf=expand("{outdir}/vcf_calls/{name}_diff.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		vcf=expand("{outdir}/vcf_calls/{name}_diff_coverage.vcf", name=config["name"], outdir=config["output_dir"])
	shell:
		"""
		cp {input.vcf} {output.vcf}
		"""

# Filtering on low quality positions
rule filter_coverage:
	input:
		low_qual=expand("{outdir}/temp/{{sample}}_low_qual",  outdir=config["output_dir"]),
		vcf=expand("{outdir}/vcf_calls/{name}_diff_coverage.vcf", name=config["name"], outdir=config["output_dir"]),
	output:	
		temp_vcf=temp(expand("{outdir}/temp/{{sample}}_temp.vcf", outdir=config["output_dir"]))
	threads:
		28
	message:
		"Filtering SNPs for all provided samples: Excluding variants where at least one isolate has low (<10X) coverage."
	shell:
		"""
		snpsift=/services/tools/snpeff/4.3r/snpEff/SnpSift.jar
		cat {input.vcf} | java -jar $snpsift intervals -x {input.low_qual} > {output.temp_vcf}
		cp {output.temp_vcf} {input.vcf}
		"""

temp_vcf_input = [''.join((output,'/temp/',file,'_temp.vcf')) for file in IDS]

# Creating a file with all the positions for the VCF
rule create_all_position_bed_file:
	input:
		temp=temp_vcf_input,
		vcf=expand("{outdir}/vcf_calls/{name}_diff_coverage.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		temp(expand("{outdir}/temp/{name}_positions.bed", name=config["name"], outdir=config["output_dir"]))
	shell:
		"""
		cat {input.vcf} | grep -v "#" | cut -f 1,2 > {output}
		"""

# Preparing AFs for each sample at each position of the VCF file
rule AF_by_position_calculation:
	input:
		temp_vcf=expand("{outdir}/temp/{{sample}}_temp.vcf", outdir=config["output_dir"]),
		bed=expand("{outdir}/temp/{name}_positions.bed", name=config["name"], outdir=config["output_dir"]),
		ref=expand("{outdir}/raw_vcf_calls/{{sample}}/reference/ref.fa", outdir=config["output_dir"]),
		bam=expand("{outdir}/raw_vcf_calls/{{sample}}/{{sample}}.bam", outdir=config["output_dir"])
	output:
		temp=expand("{outdir}/temp/{{sample}}_temp_pos", outdir=config["output_dir"]),
		temp_mpileup=temp(expand("{outdir}/temp/{{sample}}_mpileup.log", outdir=config["output_dir"]))
	params:
		depth=10000
	shell:	
		"""
		module load samtools/1.8

		(while read chr pos
		do
			x=`samtools mpileup -r ${{chr}}:${{pos}}-${{pos}} -d {params.depth} -f {input.ref} {input.bam} | awk '{{print $4"\t"$5}}' | sed 's/\.+/+/g' | sed 's/\.-/-/g' | awk 'BEGIN{{FS="\t"}} {{gsub("[^.,]","",$2); if ($1 >0) print length($2)/$1; else print 0}}'`
			echo $chr $pos $x
		done < {input.bed}  > {output.temp}) 2> {output.temp_mpileup}

		"""

input_AF = [''.join((output,'/temp/',file,'_temp_pos')) for file in IDS]

# Filtering the VCF file by AF. All variants with AF>80% are excluded
rule filter_by_AF:
	input:
		temp_pos=input_AF,
		vcf=expand("{outdir}/vcf_calls/{name}_diff_coverage.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		vcf=expand("{outdir}/output_files/{name}_final.vcf", name=config["name"], outdir=config["output_dir"]),
		positions=temp(expand("{outdir}/temp/{name}_exclude_positions", name=config["name"], outdir=config["output_dir"]))
	message:
                "Filtering SNPs for all provided samples: Excluding variants where all isolates have >80% alternative allele frequency"
	params:
		len=sample_no,
		len2=sample_no2
	shell:
		"""
		snpsift=/services/tools/snpeff/4.3r/snpEff/SnpSift.jar
		
		cat {input.temp_pos} | awk '$3<0.8' | cut -d " " -f 1-2 | sort | uniq -c | awk -v num={params.len} -v num2={params.len2} '$1==num || $1==num2 {{print $2,$3-1,$3}}'  | tr " " "\t" > {output.positions}

		cat {input.vcf} | java -jar $snpsift intervals -x {output.positions} > {output.vcf}
		"""

# Split multi-vcf file to vcf file for each sample
rule split_final_vcf:
	input:
		vcf=expand("{outdir}/output_files/{name}_final.vcf", name=config["name"], outdir=config["output_dir"])
	output:
		split_vcf=temp(expand("{outdir}/temp/{{sample}}_final_split.vcf", outdir=config["output_dir"])),
	params:
		sample="{sample}"
	message:
		"Splitting final filtered VCF file to separate files for each sample"
	shell:
		"""
		module load perl/5.24.0
		module load bcftools/1.9
		module load vcftools/0.1.16
		
		vcf-subset --exclude-ref -t SNPs -c {params.sample} {input.vcf} > {output.split_vcf}
		
		"""

# Create fasta file from filtered (final) vcf for each sample
rule create_fasta_for_vcf:
	input:
		vcf=expand("{outdir}/temp/{{sample}}_final_split.vcf", outdir=config["output_dir"]),
		ref=expand("{outdir}/raw_vcf_calls/{{sample}}/reference/ref.fa", outdir=config["output_dir"])
	output:
		fasta=expand("{outdir}/temp/{{sample}}_final_split.fasta", outdir=config["output_dir"]),
		vcf_gz=temp(expand("{outdir}/temp/{{sample}}_final_split.vcf.gz", outdir=config["output_dir"])),
		tbi=temp(expand("{outdir}/temp/{{sample}}_final_split.vcf.gz.tbi", outdir=config["output_dir"]))
	log:
		expand("{outdir}/logs/fasta_from_vcf_{{sample}}.log",outdir=config["output_dir"])
	params:
		sample="{sample}"
	message:
		"Creating a FASTA file from a VCF file"
	shell:
		"""
		module load bcftools/1.9
		
		bgzip {input.vcf}
		tabix -p vcf {output.vcf_gz}
		
		(bcftools consensus -f {input.ref} -H A {output.vcf_gz} -o {output.fasta}) 2> {log}
		sed -i 's/^>.*/>{params.sample}/' {output.fasta}
		"""	

consensus_fasta_list = [''.join((output,'/temp/',file,'_final_split.fasta')) for file in IDS]

# Create a multi-FASTA file for all samples
rule create_multi_fasta:
	input:
		fasta_list=consensus_fasta_list
	output:
		multi_fasta=expand("{outdir}/output_files/{name}_final_filtered.fa", outdir=config["output_dir"], name=config["name"]),
		temp_oneliner=temp(expand("{outdir}/output_files/{name}_final_filtered_temp_oneliner.fa", outdir=config["output_dir"], name=config["name"])),
		temp_out=temp(expand("{outdir}/output_files/{name}_final_filtered_temp.fa", outdir=config["output_dir"], name=config["name"]))
	message:
		"Creating a multi-FASTA file for all samples"
	shell:
		"""
		cat {input.fasta_list} > {output.temp_out}
	
		awk \'/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}\' {output.temp_out} | tail -n +2 > {output.temp_oneliner}		
	
		awk \'/>/  {{ id = $0 }} !/>/ {{ seq[id] = seq[id] $0 }} END  {{ for (id in seq) print id "\\n" seq[id] }}\' {output.temp_oneliner} > {output.multi_fasta}
		
		rm {input.fasta_list}
		"""


# Create snp-distance from multi-fasta file
rule calculate_snp_distances:
	input:
		multi_fasta=expand("{outdir}/output_files/{name}_final_filtered.fa", outdir=config["output_dir"], name=config["name"])
	output:
		snp_dist=expand("{outdir}/output_files/{name}_snp_dist.tsv", outdir=config["output_dir"], name=config["name"])
	params:
		snp_dists="/home/projects/cu_10047/people/misur/snp-dists/snp-dists"
	message:
		"Creating SNP distance matrix"
	shell:
		"""
		{params.snp_dists} -q {input.multi_fasta} > {output.snp_dist}
		"""

# Create RAxML tree from a multi-fasta file
rule generate_raxml_tree:
	input:
		multi_fasta=expand("{outdir}/output_files/{name}_final_filtered.fa", outdir=config["output_dir"], name=config["name"])
	output:
		raxml=expand("{outdir}/output_files/RAxML_bestTree.{name}", outdir=config["output_dir"], name=config["name"]),
		raxml_log=temp(expand("{outdir}/output_files/RAxML_log.{name}", outdir=config["output_dir"], name=config["name"])),
		raxml_result=temp(expand("{outdir}/output_files/RAxML_result.{name}", outdir=config["output_dir"], name=config["name"])),
		raxml_info=temp(expand("{outdir}/output_files/RAxML_info.{name}", outdir=config["output_dir"], name=config["name"])),
		raxml_parsimony=temp(expand("{outdir}/output_files/RAxML_parsimonyTree.{name}", outdir=config["output_dir"], name=config["name"]))
	log:
		expand("{outdir}/logs/raxml_{name}.log",outdir=config["output_dir"], name=config["name"])
	params:
		name=config["name"],
		method="GTRCAT",
		p="12345",
		dir=expand("{outdir}/output_files/", outdir=config["output_dir"])
	message:
		"Creating a phylogenetic tree based on SNPs"
	shell:
		"""
		module load raxml/8.2.11

		(raxmlHPC -n {params.name} -s {input.multi_fasta} -m {params.method} -p {params.p} -w {params.dir}) 2> {log}
		"""

if len(IDS)>0:
	first_sample_ID=IDS[0]

consensus_fasta_list = ''.join((output,'/raw_vcf_calls/',first_sample_ID,'/reference/ref.fa'))


# Create another RAxML tree with the reference genome
rule generate_raxml_tree_with_reference:
	input:
		multi_fasta=expand("{outdir}/output_files/{name}_final_filtered.fa", outdir=config["output_dir"], name=config["name"]),
		ref=consensus_fasta_list
	output:
		multi_fasta_with_ref=expand("{outdir}/output_files/{name}_final_filtered_incl_reference.fa", outdir=config["output_dir"], name=config["name"]),
		raxml=expand("{outdir}/output_files/RAxML_bestTree.{name}", outdir=config["output_dir"], name=name_with_ref),
		raxml_log=temp(expand("{outdir}/output_files/RAxML_log.{name}", outdir=config["output_dir"], name=name_with_ref)),
                raxml_result=temp(expand("{outdir}/output_files/RAxML_result.{name}", outdir=config["output_dir"], name=name_with_ref)),
                raxml_info=temp(expand("{outdir}/output_files/RAxML_info.{name}", outdir=config["output_dir"], name=name_with_ref)),
		raxml_parsimony=temp(expand("{outdir}/output_files/RAxML_parsimonyTree.{name}", outdir=config["output_dir"], name=name_with_ref)),
		reference=temp(expand("{outdir}/temp/{name}.fa", outdir=config["output_dir"], name=name_with_ref)),
		reference_oneliner=temp(expand("{outdir}/temp/{name}_oneliner.fa", outdir=config["output_dir"], name=name_with_ref))
	log:
		expand("{outdir}/logs/raxml_{name}.log",outdir=config["output_dir"], name=name_with_ref)
	params:
		name=name_with_ref,
                method="GTRCAT",
                p="12345",
                dir=expand("{outdir}/output_files/", outdir=config["output_dir"])
	message:
		"Creating a phylogenetic tree based on SNPs with the reference genome"
	shell:
		"""
		module load raxml/8.2.11

		name=`cat {input.ref} | grep ">" | head -1`
		seq=`cat {input.ref} | grep -v ">"`
		echo "$name" $'\\n' "$seq" > {output.reference}
	
		awk \'/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}\' {output.reference} | tail -n +2 > {output.reference_oneliner}

		cat {output.reference_oneliner} {input.multi_fasta} > {output.multi_fasta_with_ref}

		ref_name=`head -1 {input.ref} | sed 's/>//'`
		(raxmlHPC -o $ref_name -n {params.name} -s {output.multi_fasta_with_ref} -m {params.method} -p {params.p} -w {params.dir}) 2> {log}

		"""

# make tree produced only if there are enough samples
