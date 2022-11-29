"""
Title: SnakeMake pipeline to analyze bulk, paired-end RNA-Seq data
Author: Garth Kong

"""

import os
import sys
import glob
import pandas as pd
import plotly as plt
import plotly.graph_objects as go
from snakemake.utils import min_version
min_version("5.11")
if sys.version_info < (3, 6):
	sys.exit("Python version is less than 3.6. Your python version:", sys.version_info)

SAMPLES, = glob_wildcards("data/raw/{sample}_R1.fastq.gz")
READS = sorted(glob.glob("data/raw/*.gz"))
READS = [ os.path.basename(i).split(".")[0] for i in READS ]

def message(msg):
	sys.stderr.write("|--- " + msg + "\n")

for i in SAMPLES:
	message("Processing " + i)

def defect_mode(wildcards, attempt):
	if attempt == 1:
		return ""
	elif attempt > 1:
		return "-D"

configfile: "config.yaml"

singularity: "library://gartician/miniconda-mamba/4.12.0:sha256.7302640e37d37af02dd48c812ddf9c540a7dfdbfc6420468923943651f795591"

rule all:
	input:
		# quality control -------------------------------------------------------------------------
		expand("data/fastp/{sample}_{read}.fastq.gz", sample = SAMPLES, read = ["R1", "R2"]),
		expand("data/fastqc/{reads}_fastqc.html", reads = READS),
		expand("data/fastq_screen/{reads}_screen.txt", reads = READS),
		expand("data/preseq/estimates_{sample}.txt", sample = SAMPLES),
		expand("data/preseq/lcextrap_{sample}", sample = SAMPLES),
		# "data/multiqc/multiqc_report.html",
		# "data/fraglen.html",
		# "data/frip.html",
		# read alignment --------------------------------------------------------------------------
		expand("data/star/{sample}_bam/Aligned.sortedByCoord.out.bam", sample = SAMPLES),
		expand("data/bigwig/{sample}.bw", sample = SAMPLES),


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
	threads: 12
	shell:
		"bamCoverage -b {input[0]} -o {output} -p {threads} --normalizeUsing CPM --binSize 10 --smoothLength 50"



# DESeq2 ------------------------------------------------------------------------------------------