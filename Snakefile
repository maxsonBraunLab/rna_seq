"""
Title: SnakeMake pipeline to analyze bulk, paired-end RNA-Seq data
Author: Garth Kong

"""

import os
import sys
from glob import glob
from pandas import read_csv
from snakemake.utils import min_version
from yaml import safe_load

min_version("5.11")

if sys.version_info < (3, 6):
	sys.exit("Python version is less than 3.6. Your python version:", sys.version_info)

SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")

configfile: "config.yaml"

# singularity: "/home/groups/MaxsonLab/software/singularity-containers/4.12.0_sha256.7302640e37d37af02dd48c812ddf9c540a7dfdbfc6420468923943651f795591.sif"

def message(msg):
	sys.stderr.write("|--- " + msg + "\n")

for i in SAMPLES:
	message("Processing " + i)

def define_reads():
	READS = sorted(glob("data/raw/*.gz"))
	READS = [ os.path.basename(i).split(".")[0] for i in READS ]
	return READS

READS = define_reads()

def defect_mode(wildcards, attempt):
	if attempt == 1:
		return ""
	elif attempt > 1:
		return "-D"

def define_contrasts(file = config["CONTRASTS"]):
	contrasts = read_csv(file, sep = "\t", header = None)
	contrast1 = contrasts[0]
	contrast2 = contrasts[1]
	return contrast1, contrast2

contrast1, contrast2 = define_contrasts()

def detect_singularity():
	from sys import argv
	cmd = " ".join(sys.argv)
	singularity_flag = "--use-singularity"
	if cmd.find(singularity_flag) != -1:
		return "true"
	else:
		return "false"


rule all:
	input:
		# quality control -----------------------------------------------------
		expand("data/fastp/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = READS),
		expand("data/fastq_screen/{reads}_screen.txt", reads = READS),
		expand("data/preseq/estimates_{sample}.txt", sample = SAMPLES),
		expand("data/preseq/lcextrap_{sample}", sample = SAMPLES),
		"data/multiqc/multiqc_report.html",
		# read alignment ------------------------------------------------------
		expand([
			"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
			"data/star/{sample}_bam/ReadsPerGene.out.tab",
			"data/star/{sample}_bam/Log.final.out"
		], sample = SAMPLES),
		expand("data/bigwig/{sample}.bw", sample = SAMPLES),
		"data/counts/{}-raw-counts.txt".format(config["PROJECT_ID"]),
		"data/counts/{}-raw-filtered-counts.txt".format(config["PROJECT_ID"]),
		# deseq2 --------------------------------------------------------------
		"data/counts/{}-deseq2-norm.txt".format(config["PROJECT_ID"]),
		"data/counts/{}-log2-deseq2-norm.txt".format(config["PROJECT_ID"]),
		"data/deseq2/group",
		expand(["data/deseq2/pairwise/{c1}-vs-{c2}-all.txt",
				"data/deseq2/pairwise/{c1}-vs-{c2}-pca.pdf",
				"data/deseq2/pairwise/{{c1}}-vs-{{c2}}-{p}.txt".format(p = config["PADJ"]),
		], zip, c1 = contrast1, c2 = contrast2)

# pre-processing ----------------------------------------------------------------------------------

rule fastp:
	input:
		r1 = "data/raw/{sample}_R1.fastq.gz",
		r2 = "data/raw/{sample}_R2.fastq.gz"
	output:
		r1 = "data/fastp/{sample}_R1.fastq.gz",
		r2 = "data/fastp/{sample}_R2.fastq.gz"
	conda:
		"envs/fastp.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastp.sif")
	log:
		"data/logs/{sample}.fastp.json"
	threads: 8
	shell:
		"fastp "
		"-i {input.r1} "
		"-I {input.r2} "
		"-o {output.r1} "
		"-O {output.r2} "
		"--detect_adapter_for_pe "
		"--thread {threads} "
		"-j {log} "
		"-h /dev/null"

rule fastqc:
	input:
		"data/fastp/{read}.fastq.gz"
	output:
		"data/fastqc/{read}_fastqc.html"
	conda:
		"envs/fastqc.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastqc.sif")
	log:
		"data/logs/fastqc_{read}.log"
	threads: 4
	shell:
		"fastqc -t {threads} --outdir data/fastqc {input} > {log} 2>&1"

rule fastq_screen:
	input:
		fastq = "data/fastp/{read}.fastq.gz",
		config = config["FASTQ_SCREEN_CONFIG"]
	output:
		"data/fastq_screen/{read}_screen.txt"
	conda:
		"envs/fastq_screen.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "fastq_screen.sif")
	log:
		"data/logs/fastq_screen_{read}.txt"
	threads: 8
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen --conf {input.config} --force {input.fastq} > {log} 2>&1"

# read alignment ----------------------------------------------------------------------------------

rule STAR:
	input:
		fwd = "data/fastp/{sample}_R1.fastq.gz",
		rev = "data/fastp/{sample}_R2.fastq.gz"
	output:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
		"data/star/{sample}_bam/ReadsPerGene.out.tab",
		"data/star/{sample}_bam/Log.final.out"
	threads: 12
	params:
		gtf=config["GTF"],
		genome_index=config["STAR"]
	conda:
		"envs/star.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "star.sif")
	shell:
		"STAR "
		"--runThreadN {threads} "
		"--runMode alignReads "
		"--genomeDir {params.genome_index} "
		"--readFilesIn {input.fwd} {input.rev} "
		"--outFileNamePrefix data/star/{wildcards.sample}_bam/ "
		"--sjdbGTFfile {params.gtf} "
		"--quantMode GeneCounts "
		"--sjdbGTFtagExonParentGene gene_name "
		"--outSAMtype BAM SortedByCoordinate "
		"--readFilesCommand zcat "
		"--twopassMode Basic"

rule index:
	input:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam"
	output:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
	conda:
		"envs/samtools.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "samtools.sif")
	threads: 4
	shell:
		"samtools index -@ {threads} {input}"

rule preseq:
	input:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
	output:
		"data/preseq/estimates_{sample}.txt"
	conda:
		"envs/preseq.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "preseq.sif")
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_{sample}.log"
	shell:
		"preseq c_curve -B {resources.defect_mode} -l 1000000000 -P -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
	input:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
	output:
		"data/preseq/lcextrap_{sample}"
	conda:
		"envs/preseq.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "preseq.sif")
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_lcextrap_{sample}.log"
	shell:
		"preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -P -e 1000000000 -o {output} {input} > {log} 2>&1"

rule bigwig:
	input:
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
		"data/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
	output:
		"data/bigwig/{sample}.bw"
	conda:
		"envs/deeptools.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "deeptools.sif")
	threads: 12
	shell:
		"bamCoverage -b {input[0]} -o {output} -p {threads} --normalizeUsing CPM --binSize 10 --smoothLength 50"

rule multiqc:
	input:
		expand("data/fastp/{sample}_{reads}.fastq.gz", sample = SAMPLES, reads = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = READS),
		expand("data/fastq_screen/{reads}_screen.txt", reads = READS),
		expand("data/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
		expand("data/preseq/estimates_{sample}.txt", sample = SAMPLES),
		expand("data/preseq/lcextrap_{sample}", sample = SAMPLES)
	output:
		"data/multiqc/multiqc_report.html"
	params:
		singularity = detect_singularity()
	conda:
		"envs/multiqc.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "multiqc.sif")
	shell:
		"if [ '{params.singularity}' == 'true' ]; then export LC_ALL=C.UTF-8; export LANG=C.UTF-8; fi && "
		"multiqc data -f --ignore data/tmp -o data/multiqc 2>&1"

# counts table ------------------------------------------------------------------------------------

rule compile_counts:
	input:
		expand("data/star/{sample}_bam/ReadsPerGene.out.tab", sample = SAMPLES)
	output:
		"data/counts/{}-raw-counts.txt".format(config["PROJECT_ID"])
	params:
		samples = SAMPLES
	script:
		"scripts/compile_star_counts.py"

rule filter_counts:
	input:
		counts = "data/counts/{}-raw-counts.txt".format(config["PROJECT_ID"]),
		annotation = config["ANNOTATION"]
	output:
		"data/counts/{}-raw-filtered-counts.txt".format(config["PROJECT_ID"])
	script:
		"scripts/filter_counts.py"

# deseq2 ------------------------------------------------------------------------------------------

# output log2-transformed deseq2-normalized counts table.
rule deseq2_norm:
	input:
		counts = "data/counts/{}-raw-filtered-counts.txt".format(config["PROJECT_ID"]),
		md = config["DESEQ2_CONFIG"]
	output:
		norm_counts = "data/counts/{}-deseq2-norm.txt".format(config["PROJECT_ID"]),
		log2_norm_counts = "data/counts/{}-log2-deseq2-norm.txt".format(config["PROJECT_ID"])
	conda:
		"envs/deseq2.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "deseq2.sif")
	script:
		"scripts/deseq2-norm.R"

rule deseq2_pairwise:
	input:
		counts = "data/counts/{}-raw-filtered-counts.txt".format(config["PROJECT_ID"]),
		md = config["DESEQ2_CONFIG"],
		contrasts = config["CONTRASTS"]
	output:
		all_genes = "data/deseq2/pairwise/{c1}-vs-{c2}-all.txt",
		sig_genes = "data/deseq2/pairwise/{{c1}}-vs-{{c2}}-{p}.txt".format(p = config["PADJ"]),
		pca_plot = "data/deseq2/pairwise/{c1}-vs-{c2}-pca.pdf"
	params:
		model = config["MODEL"],
		padj = config["PADJ"],
		c1 = "{c1}",
		c2 = "{c2}",
		outdir = "data/deseq2/pairwise"
	conda:
		"envs/deseq2.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "deseq2.sif")
	threads: 4
	log:
		"data/logs/deseq2-pairwise-{c1}-vs-{c2}.log"
	script:
		"scripts/deseq2-pairwise.R"

rule deseq2_group:
	input:
		counts = "data/counts/{}-raw-filtered-counts.txt".format(config["PROJECT_ID"]),
		md = config["DESEQ2_CONFIG"]
	output:
		outdir = directory("data/deseq2/group")
	params:
		model = config["MODEL"],
		padj = config["PADJ"]
	conda:
		"envs/deseq2.yaml"
	singularity:
        os.path.join(config["SINGULARITY_IMAGE_FOLDER"], "deseq2.sif")
	log:
		"data/logs/deseq2-group.log"
	script:
		"scripts/deseq2-group.R"

