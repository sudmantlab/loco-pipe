Analyses using the full global SNP list
================

- [SNP calling
  (`snp_calling_global.smk`)](#snp-calling-snp_calling_globalsmk)
  - [snp_calling_global](#snp_calling_global)
- [Local PCA
  (`run_lostruct_global.smk`)](#local-pca-run_lostruct_globalsmk)
  - [split_beagle_global](#split_beagle_global)
  - [run_pcangsd_in_windows_global](#run_pcangsd_in_windows_global)
  - [summarize_pcangsd_for_lostruct_global](#summarize_pcangsd_for_lostruct_global)
  - [run_lostruct_global](#run_lostruct_global)
  - [plot_lostruct_mds_global](#plot_lostruct_mds_global)
  - [run_pcangsd_with_lostruct_outliers_global](#run_pcangsd_with_lostruct_outliers_global)
  - [plot_lostruct_outlier_pca_global](#plot_lostruct_outlier_pca_global)

## SNP calling (`snp_calling_global.smk`)

#### snp_calling_global

This rule uses ANGSD to identify single nucleotide polymorphisms (SNPs)
that segregate among all samples based on a probabilistic framework. It
also generates genotype likelihoods, depth counts, minor allele
frequency at these SNP positions, as well as a covariance matrix and a
distance matrix with all samples. This rule runs on each chromosome
separately.

- Key output files
  - `angsd/snp_calling_global/{chr}.snp_list`: a list of SNPs that
    contains information of SNP positions as well as their major and
    minor alleles (see <https://www.popgen.dk/angsd/index.php/Sites> for
    details).
  - `angsd/snp_calling_global/{chr}.mafs.gz`: minor allele frequency at
    each SNP position (see
    <https://www.popgen.dk/angsd/index.php/Allele_Frequencies> for
    details).
  - `angsd/snp_calling_global/{chr}.beagle.gz`: genotype likelihoods of
    each sample at each SNP position in beagle format (see
    <https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format>
    for details).
  - `angsd/snp_calling_global/{chr}.ibs.gz`: a randomly chosen base in
    each individual at each SNP position (see
    <https://www.popgen.dk/angsd/index.php/PCA_MDS> for details).
  - `angsd/snp_calling_global/{chr}.covMat`: a covariance matrix across
    all samples, which can be used for a principal component analysis
    (PCA) (see <https://www.popgen.dk/angsd/index.php/PCA_MDS> for
    details).
  - angsd/snp_calling_global/{chr}.ibsMat\`: a distance matrix across
    all samples, which can be used for a principal coordinate analysis
    (PCoA) (see <https://www.popgen.dk/angsd/index.php/PCA_MDS> for
    details).
  - `angsd/snp_calling_global/{chr}.pos.gz`: total depth count summed
    across all samples at each SNP position (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
  - `angsd/snp_calling_global/{chr}.depthGlobal`: depth histogram across
    all samples at SNP positions (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
  - `angsd/snp_calling_global/{chr}.depthSample`: depth histogram for
    each sample at SNP positions (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
- Hard-coded arguments (see
  <https://www.popgen.dk/angsd/index.php/SNP_calling> for details):
  - `bam {input.bamlist}`: path to the global bamlist created by rule
    `get_bamlist_global`.
  - `-doGlf 2`: export beagle-formatted likelihood files.
  - `-doMaf 1`: return minor allele frequency at each site
  - `-doMajorMinor 1`: infer the major and the minor alleles directly
    from genotype likelihoods.
  - `-doDepth 1`: return the depth distribution for each individual and
    all individuals combined.
  - `-doCounts 1`: return count of alleles at each site
  - `-maxDepth 100000`: Sites with more than 100000 reads are counted as
    having 100000 reads in the `.depthGlobal` and `.depthSample` files.
  - `-dumpCounts 1`: print the summed depth of all individuals at each
    site in the `.pos.gz` file
  - `-doIBS 1`: randomly print a single base from each individual at
    each position
  - `-makeMatrix 1`: print out a matrix of average distance between
    pairs of individuals
  - `-doCov 1`: print out the covariance matrix which can be used for
    PCA
  - `-minDepth $MINDP maxDepth $MAXDP`: minimum and maximum read depth a
    site can have to be considered as a SNP. These values are read from
    the `depth_filter.tsv` file generated by “get_depth_filter_global”.
  - `-remove_bads 1`: remove duplicated or failed reads.
  - `-only_proper_pairs 1`: only include paired end data that both reads
    of a pair are mapped correctly
  - `-r {wildcards.chr}`: constrain the analysis on each chromosome
    separately. This ensures that all chromosomes are analyzed in
    parallel.
- Customizable arguments (see
  <https://www.popgen.dk/angsd/index.php/SNP_calling> for details):
  - `-ref`: path to the reference genome. Users can specify their
    reference genome in the config file.
  - `-P`: number of threads this rule uses. This is defaulted to be 8
    but can be modified in the config file.
  - `-minQ`: minimum sequence quality threshold. This is defaulted to be
    20 but can be modified in the config file.
  - `-minMapQ`: minimum mapping quality threshold. This is defaulted to
    be 20 but can be modified in the config file.
  - `-GL` : the genotype likelihood model that ANGSD uses. This is
    defaulted to be 1 but can be modified in the config file. 1 for
    SAMtools, 2 for GATK, 3 for SOAPsnp, and 4 for SYK (see
    <https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Options>
    for details)
  - `-SNP_pval`: SNP p-value threshold. This is defaulted to be 1e-6 but
    can be modified in the config file.
  - `-minMaf`: minimum minor allele frequency threshold for a site to be
    considered as a SNP. This is defaulted to be 0.05 can be modified in
    the config file.
  - `-minInd`: minimum number of individuals that have to have data for
    a site to be included. This is defaulted to be 50% of the total
    sample size but can be modified in the config file.
  - `-setMinDepthInd`: minimum read depth an individual must have to be
    included in the count of individuals for `-minInd`. This is
    defaulted to be 1 and can be modified in the config file.
  - `{params.extra}`: optional arguments to be passed to the SNP calling
    step in ANGSD. These arguments could add extra input files to
    analyses (e.g. ancestral genome and individual inbreed
    coefficients). They can also be used to specify additional filters
    or additional analyses. For example, users can use the “-site” flag
    here to constrain analysis on a predefined set of SNPs. Or, users
    can use the “-rmTriallelic” flag followed by a p-value threshold to
    remove sites with more than two alleles.

## Local PCA (`run_lostruct_global.smk`)

This Snakefile has seven rules that run local PCA with `lostruct` and
visualize the result. We note that the local PCA analysis is run in two
different ways by default. First, it considers each chromosome
separately. Then, it combines all chromosomes together. The “separate”
approach is often good for identifying smaller outlier windows showing
unique PCA patterns in each chromosome, while outliers along the same
MDS axis do not necessarily share similar PCA patterns. The “combined”
approach is often good for identifying larger outlier windows
(e.g. large inversions several Mbp in size) or outlier windows located
at different windows that share similar PCA patterns (e.g. multiple
independent affected by strong divergent selection in a system with high
gene flow).

#### split_beagle_global

This rule splits chromosome-level beagle formatted genotype likelihoods
into windows with fixed number of SNPs. This rule runs on each
chromosome separately.

- Key output files
  - `{basedir}/lostruct/global/split_beagle/{chr}.w*.beagle.gz`:
    genotype likelihood in windows in beagle format
- Customizable arguments:
  - `snp_window_size`: number of SNPs in a window. This is defaulted to
    be 100 and can be modified in the config file. The smaller this
    number is, the smaller each window is, and the more windows there
    will be. Having smaller windows will result in higher computational
    burden for the downstream, but has more power in detecting small
    outlier regions.

#### run_pcangsd_in_windows_global

This rule runs PCAngsd at each window. This rule runs on each chromosome
separately.

- Key output files
  - `lostruct/global/run_pcangsd_in_windows/{chr}.w{window_id}.cov`:
    covariance matrix at each genomic window on each chromosome
- Customizable arguments:
  - `minmaf`: minimum allele frequency for a SNP to be considered in the
    local PCA analysis. It is defaulted to be 0.05 but can be modified
    in the config file.

#### summarize_pcangsd_for_lostruct_global

This rule summarizes PCA results at each window as required by the
lostruct package. This rule runs on each chromosome separately.

- Key output files
  - `lostruct/global/run_pcangsd_in_windows/{chr}.pca_summary.tsv`: PCA
    summary stats as required by lostruct
  - `lostruct/global/run_pcangsd_in_windows/{chr}.snp_position.tsv`:
    start and end position of each genomic window

#### run_lostruct_global

This rule runs the lostruct package using PCA summaries as input and
outputs a distance matrix, which describes the dissimilarity between
windows. A multidimensional scaling (MDS) is then performed on the
distance matrix.

- Key output files
  - `lostruct/global/run_lostruct/{chr}.mds.tsv`: MDS result when each
    chromosome is considered separately.
  - `lostruct/global/run_lostruct/combined.mds.tsv`: MDS result when all
    chromosomes are combined.

#### plot_lostruct_mds_global

This rule plots the top MDS axes and saves the outlier windows in text
files.

- Key output files
  - `figures/lostruct/global/separated.mds.png`: MDS plot when each
    chromosome is considered separately.
  - `figures/lostruct/global/combined.mds.png`: MDS plot when all
    chromosomes are combined.
- Customizable arguments
  - `z_cutoff`: the z score cutoff for windows to be considered outliers
    along each MDS axis. It is defaulted to be 3 but can be modified in
    the config file.

#### run_pcangsd_with_lostruct_outliers_global

This rule combines outlier windows in each direction along the same MDS
axis and performs a consensus PCA with PCAngsd.

- Key output files
  - `lostruct/global/run_pcangsd_with_lostruct_outliers/{chr}.mds_{mds_axis}.{sign}.cov`:
    covariance matrix generated from SNPS in outlier windows on each
    chromosome along each direction (plus or minus) of each MDS axis.
  - `lostruct/global/run_pcangsd_with_lostruct_outliers/combined.mds_{mds_axis}.{sign}.cov`:
    covariance matrix generated from SNPS in outlier windows across all
    chromosomes along each direction (plus or minus) of each MDS axis.

#### plot_lostruct_outlier_pca_global

This rule plots the result of the outlier windows’ consensus PCA.

- Key output files
  - `figures/lostruct/global/separated.pca.pdf`: consensus PCA plots
    when each chromosome is considered separately. Each page in the pdf
    file is one direction along a MDS axis, and if a chromosome has
    outliers along this direction, a consensus PCA at these outlier
    windows is plotted. Each point is a sample, and they are colored by
    population.
  - `figures/lostruct/global/combined.pca.png`: consensus PCA plots when
    all chromosomes are combined. Each facet contains a PCA plot at
    outlier windows in one direction along a MDS axis. Each point is a
    sample, and they are colored by population.