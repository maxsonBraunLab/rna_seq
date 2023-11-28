# rna_seq

[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)
[![snakemake minimum](https://img.shields.io/badge/snakemake->=5.32-<COLOR>.svg)](https://shields.io/)
![Maintainer](https://img.shields.io/badge/maintainer-gartician-blue)

This SnakeMake pipeline processes aligns bulk RNA-Seq datasets using STAR and finds differential gene expression between conditions using DESeq2. 

## 1. Configure your project directory

```bash
# create a new working directory and clone this repository
mkdir my_project
cd my_project
git clone https://github.com/maxsonBraunLab/rna_seq.git
cd rna_seq

# create a folder to symlink files into
mkdir -p data/raw
ln -s /absolute/path/to/files/EXP0000001_sample_R1_S001.fastq.gz
ln -s /absolute/path/to/files/EXP0000001_sample_R2_S001.fastq.gz
...

# rename symlinks to remove extraneous strings and match the following format:  {sample}_{R1|R2}.fastq.gz
mv EXP0000001_sample_R1_S001.fastq.gz sample_R1.fastq.gz
mv EXP0000001_sample_R2_S001.fastq.gz sample_R2.fastq.gz

# make scripts executable
chmod +x src/*.py src/*.sh *.sh
```

__IMPORTANT:__ Please double check your samples are in the following format before moving forward: `{sample}_{R1|R2}.fastq.gz`


## 2. Prepare your pipeline configurations

To tailor the analysis to your needs, there are a series of files to edit before running the pipeline. Example formatting for the metadata, contrasts, groups, and replicates files can be found in the `config` folder.

1. `config.yaml`
	* This file specifies pipeline runtime configurations, such as which STAR genomes to align to, feature biotypes to select (e.g. protein coding genes), and differential analysis schemes. 
2. `config/metadata.txt`
	* This is a REQUIRED 2-column table (tab-separated format, with SampleID and Condition columns) containing all unique sample IDs and the experimental condition that each sample belongs to. 
	* The sample IDs here should be the same as the sample names used in the formatted FASTQ symlinks: `{samplename}_{R1|R2}.fastq.gz`. For example, if the FASTQ files for a given sample are CSF3R-rep3_R1.fastq.gz and CSF3R-rep3_R2.fastq.gz, then the sample ID in the metadata.txt file should be CSF3R-rep3.
3. `config/contrasts.txt`
   * This is a REQUIRED 2-column table (tab-separated format, no header or column names) containing pairwise condition combinations to assess. The order of conditions matter (column1-vs-column2) to facilitate biological interpretation. 
4. `config/groups.txt`
   * An OPTIONAL 1-column table containing all the conditions used for group DESeq2 analysis. The order of conditions affect the plot order in the output heatmap. By default, use all the conditions in the metadata.
5. `config/replicates.yaml`
   * An OPTIONAL YAML file that shows how to collapse biological replicates for DESeq2 group analysis (mainly affects the heatmap). Ignored if MERGE_REPLICATES is False.


## 3. Prepare your software environment(s)

This pipeline relies on either the package manager [Conda](https://docs.conda.io/en/latest/) or the container platform [Singularity](https://apptainer.org/user-docs/master/) to create reproducible runtime environments. 


### 1. Create an environment containing Snakemake

To run the pipeline, we need a separate Conda environment with a clean install of SnakeMake - this environment can be used to invoke other maxsonBraunLab SnakeMake pipelines (including those that are containerized with Singularity), so we only need to install it once. 

If you do not have a Conda environment with Snakemake installed, please use one of the following options to create one:

* __Option A:__ Follow the instructions in the maxsonBraunLab [Snakemake setup repository](https://github.com/maxsonBraunLab/snakemake_setup) to create a Snakemake environment from a pre-made environment configuration file (recommended).
* __Option B:__ Use [Mamba](https://github.com/mamba-org/mamba) to create an environment with Snakemake installed as shown below (requires Mamba installation prior to running):

```bash
# create an env with only snakemake inside of it
mamba create -n snakemake -c bioconda -c conda-forge snakemake
```

### 2. Activate Snakemake environment

You will need to activate your Snakemake environment prior to running the pipeline. To do this, run:

```bash
conda activate <name_of_snakemake_environment>
```

### 3. Create required pipeline environments

__Note:__ If you are running the containerized Singularity version of the pipeline, then you can skip the rest of the Conda environment installation steps below and proceed to the ["Reproducible results with SnakeMake + Singularity"](#4-optional-reproducible-results-with-snakemake--singularity) section.

If you are running the non-containerized Conda version of the pipeline, then you will need to create all the environments required to run the pipeline. To avoid having to install these environments every time you run the pipeline on a new dataset, you can save all your environments to a specific location as shown below. 

```bash
# while in the working directory (where the Snakefile is), create all the environments required in the pipeline.
# specifically save them to this prefix folder so we only install everything once.
conda_prefix="${CONDA_PREFIX_1}/envs"

snakemake -j 1 --use-conda --conda-prefix $conda_prefix --conda-create-envs-only
```

## 4. (Optional) Reproducible results with SnakeMake + Singularity

To ensure the reproducibility of your results and reduce environment setup issues, we recommend running a SnakeMake workflow using Singularity containers. These containers standardize the installation of bioinformatics software (e.g. bowtie2, samtools, deseq2). 

If you would like to run the pipeline using Singularity containers instead of Conda, please follow additional setup instructions below.

### SnakeMake + Singularity Setup

If you have access to the MaxsonLab storage space on Exacloud, then you can skip the first setup step below and use the default `SINGULARITY_IMAGE_FOLDER` path specified in the `config.yml` file to access the containers for this pipeline. 

If you have previously run this pipeline with Singularity and already built the needed containers, then you can set the `SINGULARITY_IMAGE_FOLDER` path in the `config.yml` to the folder where your container images are stored.

#### 1. (Optional) Build Singularity containers

Do this step if you do **not** have access to the MaxsonLab storage space on Exacloud and have not run this pipeline with Singularity before. You will need to build the necessary containers from the definition files provided in the `singularity_definition_files` folder. 

To build containers without requiring root access on Exacloud, you will need to create a [Sylabs](https://cloud.sylabs.io/) account (you can use your Github account to log in). After logging in, navigate to `Dashboard > Access Tokens` and create a new access token. Make sure to copy and save the token into a secure place. 

This token will allow you to access the Sylabs remote builder tool from the command line. Note that every user is limited to 500 minutes of build time per month.

After generating a Sylabs access token, you will need to log into Sylabs from the Exacloud command line. To do this:

```bash
# get onto an interactive/compute node
srun -p light --time=3:00:00 --pty bash

# load the singularity module (only available on a compute node)
module load /etc/modulefiles/singularity/current

# input your access token when prompted
singularity remote login
```

Now you're ready to start building containers. To do this, navigate to the main folder of the pipeline (where the Snakefile is) and run the `singularity_build_remote.sh` script as follows:

```bash
# create folder to store build logs
mkdir -p jobs/singularity_build_remote

# make sure to follow any additional instructions in the script file before executing

# provide the path to a folder where you want to store your container images 
# (don't include a slash "/" at the end of path)
sbatch singularity_build_remote.sh <path_to_output_folder>
```

#### 2. Set Singularity cache directory

By default, Singularity will create and use a cache directory in your personal user root folder (i.e. in `/home/users/<username>`). This may create problems as there is limited space in a user's root folder on Exacloud. To avoid issues with space in your root folder, you can set the Singularity cache directory path to a folder in your lab group directory like this:

```bash
# make a cache folder inside your lab user folder 
mkdir /home/groups/MaxsonLab/<your_user_folder>/singularity_cache

# make the path to the cache folder accessible to other processes
export SINGULARITY_CACHEDIR="/home/groups/MaxsonLab/<your_user_folder>/singularity_cache"
```

If you are an experienced user, you can add the `export SINGULARITY_CACHEDIR=...` line above to your `~/.bashrc` file. Otherwise, run the `export SINGULARITY_CACHEDIR=...` command before executing the pipeline.



## 5. Set up SLURM integration (for batch jobs)

Please follow the instructions in the "Snakemake + SLURM integration" section below if you are running the pipeline as a batch job and don't yet have a [SLURM profile](https://github.com/Snakemake-Profiles/slurm) set up. The SLURM profile will configure default settings for SnakeMake to interact with SLURM. More information can be found [here](https://github.com/maxsonBraunLab/slurm).

**NOTE:** If you already have a SLURM profile set up to run Snakemake with Conda (i.e., includes settings like use-conda, conda-prefix) but would like to run Snakemake with SLURM and Singularity integration, please follow the instructions in the "Snakemake + SLURM + Singularity integration" section below.

### Snakemake + SLURM integration

Download the `slurm` folder from the maxsonBraunLab [repository](https://github.com/maxsonBraunLab/slurm) and copy the entire thing to `~/.config/snakemake`. 

Your file configuration for SLURM should be as follows:

```
~/.config/snakemake/slurm/<files>
```

Change the file permissions for the scripts in the `slurm` folder so that they are executable. To do this, run:

```bash
chmod +x ~/.config/snakemake/slurm/slurm*
```

### Snakemake + SLURM + Singularity integration

**NOTE:** The `~/.config/snakemake/slurm/config.yaml` file contains settings for SnakeMake to interact with SLURM and, optionally, Conda or Singularity. If you already have an exisiting SLURM profile configured to run Snakemake with Conda (i.e., includes settings like use-conda, conda-prefix), then you will need to create a separate profile for running Snakemake with Singularity. To do this: 

1. Copy contents of base slurm profile into another folder for slurm_singularity profile:

```bash
cp -r ~/.config/snakemake/slurm ~/.config/snakemake/slurm_singularity
```

2. Make profile scripts executable:

```bash
chmod +x ~/.config/snakemake/slurm_singularity/slurm*
```

3. Remove any conda-specific settings (e.g. use-conda, conda-prefix, etc.) from `~/.config/snakemake/slurm_singularity/config.yaml`

4. (Optional) Add the following lines to end of `~/.config/snakemake/slurm_singularity/config.yaml` file:

```
use-singularity: True
keep-going: True
rerun-incomplete: True
printshellcmds: True
```


## 6. Run the pipeline

The pipeline can either be run using Conda package management (option 1) or using Singularity containers (option 2, recommended). For either execution option, the pipeline can be run in interactive mode or in batch mode.

Interactive mode is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process.

Batch mode execution with SLURM is appropriate for running many computationally-intensive tasks such as read alignment. 

### Option 1: Snakemake + Conda Execution

A "dry-run" can be accomplished to see what and how files will be generated by using the command:

```bash
snakemake -nrp
```

To invoke the pipeline, please use either of the two options below:

```bash
# Option A: run in interactive mode. Recommended for running light jobs.
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j <n cores> --use-conda --conda-prefix $CONDA_PREFIX_1/envs

# Option B: run in batch mode. Recommended for running intensive jobs.
sbatch run_pipeline_conda.sh
```

For users running the pipeline in batch mode, `run_pipeline_conda.sh` is a wrapper script that contains the following command:

```bash
snakemake -j <n jobs> --use-conda --conda-prefix $CONDA_PREFIX_1/envs --profile slurm --cluster-config cluster.yaml
```

Additional setup instructions are provided in the wrapper script.

You can standardize further arguments for running the pipeline in batch mode using the following [instructions](https://github.com/Snakemake-Profiles/slurm). The maxsonBraunLab repository [slurm](https://github.com/maxsonBraunLab/slurm) contains further instructions to set up a SnakeMake slurm profile.

<!--
#### Local execution

You can run the pipeline using an interactive node like this:

```bash
srun --cores=20 --mem=64G --time=24:00:00 --pty bash
conda activate snakemake
snakemake -j <n cores> --use-conda
```

This is sufficient for small jobs or running small parts of the pipeline, but not appropriate for the entire process.

#### SLURM execution

To run the pipeline using batch mode use the following command:

```bash
snakemake -j 64 \
	--use-conda \
	--conda-prefix $conda_prefix \
	--cluster-config cluster.yaml \
	--profile slurm 
```

This will allow SnakeMake to submit up to 64 jobs at once. The conda environments will be drawn from the $conda_prefix folder.  SLURM execution is appropriate for running many computationally-intensive programs (e.g. read alignment).
-->


### Option 2: Snakemake + Singularity Execution

More Singularity documentation on Exacloud can be found [here](https://wiki.ohsu.edu/display/ACC/Exacloud%3A+Singularity). 

A "dry-run" can be accomplished to see what and how files will be generated by using the command:

```bash
snakemake -nrp
```

To invoke the pipeline, please use either of the two options below. 

**NOTE:** Make sure to use double quotes for the `--bind` argument, and insert an integer for the `-j` flag. The `--bind` argument binds the host (Exacloud) paths to the container to access the genome indices and the path to the raw sequencing files. When Snakemake is executed directly on an interactive terminal, the `-j` flag represents the max number of cores to use. When executed via a batch script, the `-j` flag represents the max number of jobs to run at a time. 

**Option A: Singularity + interactive run**

```bash
# get onto an interactive/compute node
srun --time=24:00:00 --pty bash

# activate your environment with snakemake
conda activate <snakemake_env>

# load the singularity module
module load /etc/modulefiles/singularity/current

# set folder paths
# fastq_folder should be absolute/full path to folder containing raw FASTQ files (not the symlinks)
indices_folder="/home/groups/MaxsonLab/indices"
fastq_folder="/set/full/path/here"

# run pipeline
snakemake -j <n_cores> \
--use-singularity \
--singularity-args "--bind $indices_folder,$fastq_folder"
```

**Option B: Singularity + slurm (batch) run**

For users running the Singularity version of the pipeline in batch mode, `run_pipeline_singularity.sh` is a wrapper script for the pipeline. You will need to add the appropriate FASTQ folder path to the script prior to running. Additional instructions are provided in the wrapper script.

```bash
# run pipeline
sbatch run_pipeline_singularity.sh
```


<!--
To utilize Singularity in your analysis, log in to an interactive node and load the Singularity module first like this:

```bash
# request an interactive node
srun -p light --time=36:00:00 --pty bash

# re-activate your environment with snakemake
conda activate <snakemake-env>

# load the singularity program
module load /etc/modulefiles/singularity/current
```

More Singularity documentation on Exacloud can be found [here](https://wiki.ohsu.edu/display/ACC/Exacloud%3A+Singularity). If it is your first time running the pipeline, and especially when  using Singularity, we must install all the conda environments using the  following command:

```bash
indices_folder="/home/groups/MaxsonLab/indices"
conda_folder="${CONDA_PREFIX_1}/envs"
fastq_folder="/home/groups/MaxsonLab/input-data2/path/to/FASTQ/folder/"

snakemake -j 1 \
	--verbose \
	--use-conda \
	--conda-prefix $conda_folder \
	--use-singularity \
	--singularity-args "--bind $indices_folder,$conda_folder,$fastq_folder" \
	--conda-create-envs-only
```

The above code snippet will take about an hour or more to  set up, but is a one-time installation. After creating the conda  environments and configuring the pipeline, we can invoke the pipeline in the same shell like this:

```bash
# Singularity + interactive run
snakemake -j <n cores> \
	--use-conda \
	--conda-prefix $conda_folder \
	--use-singularity \
	--singularity-args "--bind $indices_folder,$conda_folder,$fastq_folder"

# Singularity + slurm run
snakemake -j <n jobs> \
	--use-conda \
	--conda-prefix $conda_folder \
	--use-singularity \
	--singularity-args "--bind $indices_folder,$conda_folder,$fastq_folder" \
	--profile slurm \
	--cluster-config cluster.yaml
```

NOTE: make sure to use double quotes and insert an integer for the -j flag.

The above command will install the pipeline's conda environments into the `conda-prefix` directory - this means that conda environments are actually not stored INSIDE the container. The `--bind` argument binds the host (Exacloud) paths to the container to access the genome indices, conda prefix, and the path to the raw sequencing files. The `--profile slurm` will configure default settings for SnakeMake to interact with SLURM - more information can be found [here](https://github.com/maxsonBraunLab/slurm). Feel free to create another [snakemake profile](https://wiki.ohsu.edu/display/ACC/Exacloud%3A+Singularity) that has its own set of singularity arguments for added convenience.

-->


# Output

```
data
├── bigwig			-- genome tracks for each sample
├── counts			-- raw and deseq2-normalized counts tables
├── deseq2			-- group and pairwise analysis of DESeq2 results
├── fastp			-- read trimming results from fastp
├── fastqc			-- FASTQ quality for each sample
├── fastq_screen	-- contamination screen for each sample
├── logs			-- log files for each rule split by sample
├── multiqc			-- multiqc report summarizing QC metrics each sample
├── preseq			-- library complexity for each sample
├── raw				-- raw sequencing reads (symlinks highly recommended)
├── star			-- BAM files for each sample
└── tmp				-- tmp folder for various programs to minimize writing to /tmp in cluster nodes.
```

# Generic Pipeline Structure

```
.
├── cluster.yaml	-- configuration file for SLURM resources
├── config			-- directory for further configuration files
├── config.yaml		-- main config file for the pipeline
├── data			-- contains raw and processed data
├── envs			-- conda environments in YAML format
├── jobs			-- cluster jobs belong here, great for debugging
├── scripts			-- series of scripts used in the pipeline
└── Snakefile		-- core of this SnakeMake pipeline
```

# Methods

![](config/rulegraph.svg)
