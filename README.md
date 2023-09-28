Loco-pipe: population-level analysis of low-coverage whole genome
sequencing data
================

Here, we present you Loco-pipe, a streamlined pipeline that could
perform exploratory analyses on high-throughput genomic data. This
pipeline will be useful to researchers who want to explore the
population structure within the total samples. We include software (etc
Ohana, PCAngsd) to calculate Fst, heterozygosity and population
admixture. Some advantages of this pipeline include: easy management by
modifying a single config file, detailed explanations embedded within
all functions, and feasibility of running the pipeline with a single
line of code. Because Loco-pipe is written in Snakemake, a workflow
management system, this pipeline also automatically inherits benefits of
Snakemake, such as ability to continue from the last failed job and to
be scalable without modifying the pipelines.

### Pipeline flowchart

![](overall_pipeline.png)

### Things you need to do before running loco-pipe on your data

1.  Build a right data structure: You **MUST** have a folder named
    exactly as “docs” within the base directory. In that directory,
    Loco-pipe will save future outputted files.

2.  Provide needed files: Within the “docs” folder you are required to
    provide two things for the pipeline to run: a **metadata table** and
    a **chromosome table**. Both need to be tsv files.

- Metadata table: This table should contain a minimum of three columns.
  A column named exactly as “sample_name”; a column named exactly as
  “bam” which stores the full pathway of all your bam files, and a
  column that stores the population-level traits you wish to segregate
  the specimens with. You will entry the name of the third column into
  the config file.

- Chromosome table: This table should contain two columns. The first
  column records the names of all chromosomes/scaffolds/contigs (the
  same names appear in the reference genome), and the second column
  records shortened names which you want to show on the Fst plot.

3.  To create the loco-pipe environment

If you do not have Conda downloaded on your computer, please follow the
[Conda
website](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
to do so. Then use the following code to download the loco-pipe
environment.

``` bash
## Download Mambdaforge. Mambaforge allows you to download with a much faster speed.
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

bash Mambaforge-$(uname)-$(uname -m).sh
## Then click "proceed" to continue.
## Once the Mambaforge is downloaded, create the loco-pipe environment with the given config file.

mamba env create -f /address/to/your/config/file/locopipe.yaml
```

4.  To download PCAngsd

We encourage you to have a separate folder for PCAngsd and follow the
steps below to finish installing. (Mambaforge is needed at this step; if
you have not downloaded mambaforge yet, please refer to the instruction
above.)

``` bash
cd /the/pathway/to/your/PCAngsd/folder  

# Download the software from the PCAngsd github website.
git clone https://github.com/Rosemeis/pcangsd.git  
cd pcangsd  
# create an environment for PCAngsd 
mamba env create -f environment.yml  
conda activate pcangsd  
python setup.py build_ext --inplace  
pip3 install -e  
conda deactivate  
```

5.  Config files

- loco-pipe config: This config file stores information that could be
  read by the loco-pipe. You should change the elements under the
  “global” tab based on your own research data and modify the
  hyper-parameters below other tabs depending on your computing power.
  This file is also designed to allow you control what features you want
  to include in your own runs. For example, if you wish not to have
  calculate Fst, you can edit the ‘get_Fst: true’ to ‘get_Fst: true’,
  and the “get_Fst” switch will be turned off.
- slurm config: This config could divide the computing units for
  different rules to improve running.
