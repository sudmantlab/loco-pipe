Pipeline preparation and depth filter determination
================

- [Pipeline preparation
  (`pipeline_prep.smk`)](#pipeline-preparation-pipeline_prepsmk)
  - [get_bamlist_global](#get_bamlist_global)
  - [get_bamlist_local](#get_bamlist_local)
  - [get_chr_list](#get_chr_list)
  - [index_bam](#index_bam)
- [Sequencing depth calculation
  (`get_depth_global.smk`)](#sequencing-depth-calculation-get_depth_globalsmk)
  - [get_depth_global](#get_depth_global)
- [Depth filter determination
  (`get_depth_filter_global.smk`)](#depth-filter-determination-get_depth_filter_globalsmk)
  - [get_depth_filter_global](#get_depth_filter_global)

## Pipeline preparation (`pipeline_prep.smk`)

This Snakefile contains four rules that prepare key input files for the
pipeline.

#### get_bamlist_global

This rule generates a text file that contains the full paths to all bam
files.

- Key output files:
  - `docs/bamlist.txt`: a text file containing the full paths to all bam
    files.

#### get_bamlist_local

This rule generates a text file that contains the full paths to all bam
files for each population (or other grouping variables that users
choose, i.e. “local”)

- Key output files:
  - `docs/{population}_bamlist.txt`: a text file containing the full
    paths to all bam files for each “population”.

#### get_chr_list

This rule extracts the first column of the chromosome table to form a
chromosome list. The resulting chromosome list is used to restrict
analyses to a predetermined set of chromosomes for certain analyses
(e.g. heterozygosity estimation).

- Key output files:
  - `docs/chr_list.txt`: a text file containing the names of chromosomes
    included in the chromosome table that users supplied.

#### index_bam

This rule uses samtools to index all the input bam files.

- Key output files:
  - `{bam}.bai`: an index file of each bam file

## Sequencing depth calculation (`get_depth_global.smk`)

#### get_depth_global

This rule uses ANGSD to count the read depth at every site summed across
all samples. This rule runs on each chromosome separately.

- Key output files (see
  <https://www.popgen.dk/angsd/index.php/Allele_Counts> for detailed
  examples):
  - `angsd/get_depth_global/{chr}.pos.gz`: total depth count summed
    across all samples at each position.
  - `angsd/get_depth_global/{chr}.depthGlobal`: depth histogram across
    all samples.
  - `angsd/get_depth_global/{chr}.depthSample`: depth histogram for each
    sample.
- Hard-coded arguments (see
  <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details):
  - `-doCounts 1`: return count of alleles at each site
  - `-dumpCounts 1`: print the summed depth of all individuals at each
    site in the `.pos.gz` file
  - `-doDepth 1`: return the depth distribution for each individual and
    all individuals combined.
  - `-maxDepth 100000`: Sites with more than 100000 reads are counted as
    having 100000 reads in the `.depthGlobal` and `.depthSample` files.
  - `-r {wildcards.chr}`: constrain the analysis on each chromosome
    separately. This ensures that all chromosomes are analyzed in
    parallel.
  - `-remove_bads 1`: remove duplicated or failed reads.
  - `-only_proper_pairs 1`: only include paired end data that both reads
    of a pair are mapped correctly
- Customizable arguments:
  - `-ref`: path to the reference genome. Users can specify their
    reference genome in the config file.
  - `-P`: number of threads this rule uses. It is defaulted to be 8 but
    can be modified in the config file.
  - `-minQ`: minimum sequence quality threshold. It is defaulted to be
    20 but can be modified in the config file.
  - `-minMapQ`: minimum mapping quality threshold. It is defaulted to be
    20 but can be modified in the config file.

## Depth filter determination (`get_depth_filter_global.smk`)

#### get_depth_filter_global

This rule reads in the depth histograms of all chromosomes
(i.e. `depthGlobal` files generated by `get_depth_global`), fits them
into a normal distribution, and establishes min and max filters based on
the mean and standard deviation of the fitted distribution. It also
outputs a density plot along with a tsv file on read depth distribution
and min and max filters.

- Key output files:
  - `angsd/get_depth_global/depth_filter.tsv`: a text file containing
    the minimum and maximum depth filters to be used in further
    analyses, as well as the mean and standard deviation of the fitted
    distribution.
  - `figures/depth/depth_filter.png`: the density plot of the read depth
    distribution. The chosen depth filters are shown as vertical lines.
    Customizable arguments:
  - `n_sd`: the maximum distance in standard deviations away from the
    mean of the fitted distribution for a site be included in the
    analyses. The higher `n_sd` is, the more relaxed depth filters will
    be, and the more sites will be retained for downstream analyses. It
    is defaulted to be 2 but can be modified in the config file.