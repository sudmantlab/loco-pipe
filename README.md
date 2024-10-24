loco-pipe: a Snakemake pipeline for low-coverage whole-genome sequencing
================

- [Key features](#key-features)
- [Currently supported
  functionalities](#currently-supported-functionalities)
- [Complete pipeline flowchart](#complete-pipeline-flowchart)
- [Before you start](#before-you-start)
- [Setting up the pipeline](#setting-up-the-pipeline)
- [Preparing the project directory and required input
  files](#preparing-the-project-directory-and-required-input-files)
- [Launching the pipeline](#launching-the-pipeline)
- [Future directions](#future-directions)
- [Citation](#citation)

**loco-pipe** is an automated Snakemake pipeline that streamlines a set
of essential population genomic analyses for **lo**w-**co**verage whole
genome sequencing (lcWGS) data.

## Key features

- Streamlining of several essential population genomic analyses
- Can be launched with a single line of code
- Incorporation of key filtering steps and best practices for
  low-coverage data
- Key results are plotted automatically for visual inspection
- Easy customization through a configuration file
- [A quick start guide with an example
  dataset](https://github.com/sudmantlab/loco-pipe/blob/main/toyfish.md)
- Extensive in-line annotation along with a detailed [user’s
  manual](manual/README.md)
- Flexible architecture that allows for the addition of new features
- Inheritance of the many benefits offered by Snakemake, including
  - High computational efficiency achieved through massive
    parallelization
  - Seamless integration with common job schedulers on computer clusters
  - The ability to automatically continue from the last failed or
    interrupted job
  - Built-in software management system and robust file structure

## Currently supported functionalities

- Depth counting
- SNP calling
- Allele frequency estimation
- Site frequency spectrum (SFS) estimation
- Fixation index (Fst) estimation
- Principal component analysis (PCA)
- Admixture analysis
- Theta estimation
- Neutrality test statistics
- Heterozygosity estimation
- Local PCA analysis

![](simplified_flowchart.png) A simplified flowchart of loco-pipe
highlighting its key functionalities. The box in dotted lines represents
user-provided input files, and boxes in solid lines represent key
analytical steps in the pipeline. Plots are generated using our example
dataset. Please see our [quick start
guide](https://github.com/sudmantlab/loco-pipe/blob/main/toyfish.md) for
detailed descriptions of the plots.

## Complete pipeline flowchart

![](complete_flowchart.png) Each box represents a Snakemake rule and is
colored based on the major groups of analyses in the form of separate
Snakefiles (shown at the top left corner). The differently colored
shades indicate the types of SNPs or sites on which the analyses are
conducted. Dashed arrows indicate key modules that can be turned on or
off in the configuration file. The four orange boxes at the top right
corner (part of the “pipeline_prep.smk” Snakefile) are the starting
points of the pipeline, and the “all” box at the bottom is the end
point. Please see our [user’s manual](manual/README.md) for detailed
descriptions of each Snakemake rule.

## Before you start

#### Reference genome

This pipeline requires a moderately contiguous reference genome for your
study system. Currently, it does not support highly fragmented reference
genomes, since most analyses are parallelized by scaffolds. Having too
many scaffolds or unscaffolded contigs will create too many parallel
jobs for Snakemake and a job scheduler to handle. We recommend that 90%
of the genome should be consisted of no more than 100 scaffolds
(i.e. L90 \< 100). Small scaffolds and contigs should be excluded from
the analysis (we will ask you to provide a list of scaffolds that you
would like to include).

#### Sequence alignment files

In addition, we assume that properly mapped and filtered bam files are
ready to be used as input files for loco-pipe. You can choose your
favorite software and/or pipeline to go from fastq to bam, but one
pipeline that we particularly recommend is
[grenepipe](https://github.com/moiexpositoalonsolab/grenepipe).

#### grenepipe

grenepipe is a Snakemake pipeline for variant calling from raw sequence
data. Although it is developed for high-coverage data, you can skip the
variant calling step by using the `all_qc`
[shortcut](https://github.com/moiexpositoalonsolab/grenepipe/wiki/Advanced-Usage#running-only-parts-of-the-pipeline)
and turning the `bcftools-stats` switch in the config file to `false`.
This way, grenepipe will stop after generating the final bam files and
their associated quality reports. (We also recommend turning the
`clip-read-overlaps` switch to `true` and setting
`VALIDATION_STRINGENCY=SILENT` for `picard MarkDuplicates` in the
configuration file.) grenepipe is very thoughtfully built and
[extensively
documented](https://github.com/moiexpositoalonsolab/grenepipe/wiki), and
it is a major inspiration for loco-pipe. Familiarizing yourself with
grenepipe will also make loco-pipe much easier to learn.

#### Quick start guide with an example dataset

If you don’t yet have your data ready for loco-pipe, and even if you do,
we highly recommend you to first follow our [quick start
guide](https://github.com/sudmantlab/loco-pipe/blob/main/toyfish.md)
which includes a heavily subsetted example dataset. It only takes
loco-pipe a few minutes to analyse the example dataset on a computer
cluster, making it much easier to learn and troubleshoot.

#### User’s manual

We also provide an extensive [user’s manual](manual/README.md) with
detailed description of each step of the pipeline for easy reference.
Tips and suggestions are also included in this document.

We **strongly** advise users to take advantage of these resources (as
well as other resources such as the [ANGSD
manual](https://www.popgen.dk/angsd/index.php/ANGSD) and our [beginner’s
guide to lcWGS](https://doi.org/10.1111/mec.16077)) to gain a good
understanding of loco-pipe, the software it uses, and lcWGS data
analysis in general while you are using loco-pipe. **Using this pipeline
as a black box can lead to spurious results and erroneous conclusions.**

## Setting up the pipeline

1.  Install
    [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
    or conda
    (<https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>)
    if you have not already done so. A fresh install of mamba with
    miniforge (<https://github.com/conda-forge/miniforge>) is highly
    recommended because mamba is much faster than conda. It is ok if you
    prefer to use conda though; just replace all occurrences of `mamba`
    with `conda` in the code below.

2.  Download `loco-pipe` from GitHub (e.g. using
    `git clone https://github.com/sudmantlab/loco-pipe.git`). We
    recommend you to download it to a folder where you store your
    software programs. We will refer to the full path of the directory
    that contains the `loco-pipe` folder as `SOFTWARE_DIR`.

3.  Create the `loco-pipe` conda environment using mamba by running
    `mamba env create -f $SOFTWARE_DIR/loco-pipe/workflow/envs/loco-pipe.yaml`
    (replace \$SOFTWARE_DIR with a real path).

4.  (Optional) If you would like to run PCA with the software
    [PCAngsd](https://github.com/Rosemeis/pcangsd) using loco-pipe, you
    **must** install PCAngsd manually as it is not yet available on
    conda. Please install it to a conda environment named
    `pcangsd_lcpipe` using the script below. Even if you already have
    PCAngsd installed on your machine, you will need to run the
    following code to ensure that the version is compatible.

    ``` bash
    # first set your working directory to a folder where you store your software programs
    cd $SOFTWARE_DIR # replace $SOFTWARE_DIR with a real path
    # download PCAngsd from Github
    git clone https://github.com/Rosemeis/pcangsd.git
    cd pcangsd
    # check out the version the loco-pipe is based on
    git checkout f90d41f7b9b245481781ae319c4a174376e5f471
    # create an environment for PCAngsd 
    mamba env create -f $SOFTWARE_DIR/loco-pipe/workflow/envs/pcangsd.yaml
    # activate the conda environment
    conda activate pcangsd_lcpipe
    # build PCAngsd
    python setup.py build_ext --inplace  
    pip3 install -e . ## if you run into issues pertaining to ssl certificates, try "pip3 install --trusted-host pypi.org -e ." instead
    # deactivate the conda environment
    conda deactivate  
    ```

5.  (Optional) If you would like to run local PCA with the
    [lostruct](https://github.com/petrelharp/local_pca) package in R
    using loco-pipe, you **must** install lostruct (in addition to
    PCAngsd, see above) manually as it is not yet available on conda.
    Please install it to a conda environment named `lostruct_lcpipe`
    using the script below.

    ``` bash
    # create a conda environment named lostruct_lcpipe and install R and some key R packages
    mamba env create -f $SOFTWARE_DIR/loco-pipe/workflow/envs/lostruct.yaml
    # activate the lostruct conda environment
    conda activate lostruct_lcpipe
    # launch R
    R
    # install lostruct
    devtools::install_github("petrelharp/local_pca/lostruct", ref = "93ad59309151e44d2d3d0d7748cdc92f6121f564")
    # quit R
    q()
    # deactivate the conda environment
    conda deactivate  
    ```

    Note: depending on your system, you may need to ensure that lostruct
    is properly installed to the `lostruct_lcpipe` environment with
    something like the following

    ``` r
    withr::with_libpaths(new = "/path/to/conda/envs/lostruct_lcpipe/lib/R/library",
                         devtools::install_github("petrelharp/local_pca/lostruct"))
    ```

## Preparing the project directory and required input files

1.  Set up the file structure.

    - First, create a base directory for your project. This folder
      should be separate from the `loco-pipe` folder. You can name it
      however you want (just avoid special characters other than dashes
      and underscores), but we will refer to the full path of this
      folder as `BASEDIR`.

    - Within `BASEDIR`, create two new folders: `docs`, and `config`.

    - You can also have your sequencing data (e.g. fastq and bam files)
      and the reference genome in separate folders in `BASEDIR`, but
      this is not required.

2.  Prepare a **sample table** and a **chromosome table**. Both of these
    should be tab separated text files. Store them in the `docs` folder
    in `BASEDIR`. You can name them however you want.

    - Sample table: This table should contain a minimum of three columns
      in no particularly order. One column should be named exactly as
      `sample_name` and it should contain sample IDs of all the samples
      to be included in the analysis. Another column should be named
      exactly as `bam` and it should store the full paths of the bam
      file for each sample. A third column should specify the grouping
      information you wish to segregate the samples by. These could be
      species, subspecies, ecotypes, populations, sampling sites, sex,
      etc. You will need to enter the name of the third column into the
      pipeline configuration file `config.yaml`. Below is an example of
      a sample table.

      | sample_name | bam                                             | species   | population  |
      |:------------|:------------------------------------------------|:----------|:------------|
      | ABLG11920-1 | /path/to/loco-pipe/toyfish/bams/ABLG11920-1.bam | sunset    | sunset      |
      | ABLG12067-1 | /path/to/loco-pipe/toyfish/bams/ABLG12067-1.bam | sunset    | sunset      |
      | ABLG11918-1 | /path/to/loco-pipe/toyfish/bams/ABLG11918-1.bam | vermilion | vermilion_1 |
      | ABLG11913-1 | /path/to/loco-pipe/toyfish/bams/ABLG11913-1.bam | vermilion | vermilion_1 |
      | ABLG9871-1  | /path/to/loco-pipe/toyfish/bams/ABLG9871-1.bam  | vermilion | vermilion_2 |
      | ABLG11795-1 | /path/to/loco-pipe/toyfish/bams/ABLG11795-1.bam | vermilion | vermilion_2 |
      | ABLG11937-1 | /path/to/loco-pipe/toyfish/bams/ABLG11937-1.bam | vermilion | vermilion_3 |
      | ABLG11940-1 | /path/to/loco-pipe/toyfish/bams/ABLG11940-1.bam | vermilion | vermilion_3 |

      > A few tips about the sample table:
      > - Sample names have to be unique, and each sample should
      >   correspond to a single bam file (i.e. we assume that if you
      >   have multiple bam files for a single sample, they have already
      >   been merged). In cases where you would like to keep different
      >   bam files from the same sample separate, e.g. for batch effect
      >   control, you will need to distinguish them by assigning them
      >   different sample names, e.g. by adding a suffix to their
      >   original names. If you do this, please also note that
      >   including multiple copies of the same sample can severely bias
      >   certain analyses, such as PCA, so we recommend you to only do
      >   this in the first iteration of the pipeline, and once the risk
      >   of batch effect is eliminated, you should merge or remove the
      >   duplicated samples.
      > - You can have multiple grouping variables in the sample table,
      >   although you may need to launch the pipeline more than once in
      >   order to conduct analyses based on different grouping
      >   variables.
      > - In practice, you may need to have multiple versions of the
      >   sample table. For example, you may include all samples in the
      >   first run, exclude some problematic one and/or change the
      >   grouping of others based on the result of the first run, and
      >   launch the pipeline again with an updated sample table.
      > - In case that you don’t have any grouping information a priori,
      >   you still need to have a fake grouping column with all samples
      >   having the same entry for plotting purposes (e.g. PCA and
      >   admixture). In this case, you will also need to turn off all
      >   population-level analyses in the configuration file (see
      >   below).

    - Chromosome table: This table should contain one or two columns.
      These columns should be **unnamed**. The first column is required,
      and it records the names of all chromosomes/scaffolds/contigs that
      you would like to include in the analysis. These should exactly
      match the names in the reference genome. The second column is
      optional, and it records shortened or alternative names of the
      chromosomes/scaffolds/contigs which you would like to show on the
      plots. If the second column is empty, the original names will be
      shown. Below is an example of a chromosome table.

      |                                             |     |
      |:--------------------------------------------|----:|
      | Sebastes_miniatus.Sebrube.F.HiC_scaffold_15 |  15 |
      | Sebastes_miniatus.Sebrube.F.HiC_scaffold_16 |  16 |

      > A few tips about the sample table:
      > - As mentioned earlier, we recommend against including too many
      >   chromosomes/scaffolds/contigs in the chromosome table (e.g. \<
      >   100).
      > - If you would like to use the plotting modules in loco-pipe to
      >   generate Manhattan-like plots automatically, we would also
      >   recommend to only include reference sequences of similar sizes
      >   (e.g. you shouldn’t include unincorporated contigs in the
      >   chromosome table when you have a chromosomal level assembly).

3.  Edit the configuration files.

    - The [pipeline config
      file](https://github.com/sudmantlab/loco-pipe/blob/main/config.yaml)
      `config.yaml`: This config file serves as an easily accessible way
      to control the behavior of loco-pipe. You can use it to 1) specify
      the location of input files, 2) include or exclude certain
      analyses from a loco-pipe run, 3) adjust the parameter settings
      for different analyses. This config file is thoroughly annotated,
      so please read it through and make edits when needed.
      **Importantly**, please copy `config.yaml` to the `config` folder
      in your project directory
      (i.e. `cp $SOFTWARE_DIR/loco-pipe/config.yaml $BASEDIR/config/config.yaml`)
      and make your edits there rather in the `loco-pipe` directory.

    - (Optional) The [cluster config
      file](https://github.com/sudmantlab/loco-pipe/blob/main/workflow/profiles/slurm/cluster_config.yaml)
      `cluster_config.yaml` under `workflow/profiles`: This config files
      specifies the resources each step of the pipeline can use on a
      computer cluster. We have extensively tested loco-pipe with a
      slurm job scheduler and an example slurm profile is included in
      loco-pipe (`workflow/profiles/slurm`), but other job schedulers
      should also work with Snakemake given an appropriate profile setup
      (see [Snakemake
      manual](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
      for details).

## Launching the pipeline

1.  Activate the `loco-pipe` environment with
    `conda activate loco-pipe`.

2.  You are now ready to launch the pipeline. We recommend to always
    start with a dry run using the `-n` flag before actually running it.

    ``` bash
    snakemake \
    --use-conda \
    --conda-frontend mamba \
    --directory $BASEDIR \
    --rerun-triggers mtime \
    --scheduler greedy \
    --snakefile $SOFTWARE_DIR/loco-pipe/workflow/pipelines/loco-pipe.smk
    ```

    When running loco-pipe on a computer cluster, make sure to use the
    `--profile` flag to specify your cluster profile. Other flags that
    may be useful are `--conda-prefix`, `--default-resources`,
    `--printshellcmds`, `--cores`, etc.

## Future directions

We plan to continue to maintain and develop loco-pipe, by incorporating
additional analyses (e.g. GWAS, dxy, LD estimation and pruning) into
this pipeline and also enabling more functionalities for the existing
software programs (e.g. ANGSD, Ohana). For certain existing analyses, we
also hope to provide more software options for users to pick from
(e.g. [winSFS](https://academic.oup.com/genetics/article/222/4/iyac148/6730749?login=false)
for SFS estimation,
[ngsAdmix](https://academic.oup.com/genetics/article/195/3/693/5935455)
for admixture analysis). Advanced SNP filtering procedures, such as the
one described in [Dallaire et
al. (2023)](https://academic.oup.com/gbe/article/15/12/evad229/7470724),
can be incorporated into loco-pipe manually (see [User’s
Manual](manual/README.md) for instructions), but may be fully automated
in the future. In choosing depth filters, we may adopt a mixed model
instead of a truncated normal distribution for curve fitting in the
future. Docker/singularity support will also be considered.

All kinds of feedback, such as bug reports and feature requests, are
welcome on the [Issues](https://github.com/sudmantlab/loco-pipe/issues)
page. We also encourage users to build on the existing infrastructure
and add more functionalities to loco-pipe in the form of [pull
requests](https://github.com/sudmantlab/loco-pipe/pulls).

## Citation

loco-pipe is published in Bioinformatics Advances:
<https://doi.org/10.1093/bioadv/vbae098> Please cite this paper if you
use loco-pipe in your research.
