Analyses using the full local SNP list
================

- [Allele frequency estimation
  (`get_maf.smk`)](#allele-frequency-estimation-get_mafsmk)
  - [get_maf](#get_maf)
- [Fst estimation (`get_fst.smk`)](#fst-estimation-get_fstsmk)
  - [get_fst](#get_fst)
  - [plot_fst](#plot_fst)

## Allele frequency estimation (`get_maf.smk`)

#### get_maf

This rule takes in population-specific bamlists and then uses ANGSD to
estimate allele frequencies at each position specified in the global SNP
list with optional population-level missingness filters (hence “local
SNP list”). It will also output population-specific genotype likelihoods
in beagle format, site allele frequency likelihoods (saf), covariance
and distances matrices, and depth counts. By supplying the global SNP
list to this rule along with the “-doMajorMinor 3” option, we ensure
that the same alleles are considered as major and minor alleles at the
same SNP position across all populations.

- Key output files
  - `angsd/get_maf/{population}.{chr}.mafs.gz`: minor allele frequency
    at each SNP position (see
    <https://www.popgen.dk/angsd/index.php/Allele_Frequencies> for
    details).
  - `angsd/get_maf/{population}.{chr}.saf.gz`: sample allele frequency
    likelihoods at each SNP position (see for
    <https://www.popgen.dk/angsd/index.php/SFS_Estimation> details)
  - `angsd/get_maf/{population}.{chr}.beagle.gz`: genotype likelihoods
    of each sample at each SNP position in beagle format (see
    <https://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#Beagle_format>
    for details).
  - `angsd/get_maf/{population}.{chr}.ibs.gz`: a randomly chosen base in
    each individual at each SNP position (see
    <https://www.popgen.dk/angsd/index.php/PCA_MDS> for details).
  - `angsd/get_maf/{population}.{chr}.covMat`: a covariance matrix
    across all samples, which can be used for a principal component
    analysis (PCA) (see <https://www.popgen.dk/angsd/index.php/PCA_MDS>
    for details).
  - `angsd/get_maf/{population}.{chr}.ibsMat`: a distance matrix across
    all samples, which can be used for a principal coordinate analysis
    (PCoA) (see <https://www.popgen.dk/angsd/index.php/PCA_MDS> for
    details).
  - `angsd/get_maf/{population}.{chr}.pos.gz`: total depth count summed
    across all samples at each SNP position (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
  - `angsd/get_maf/{population}.{chr}.depthGlobal`: depth histogram
    across all samples at SNP positions (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
  - `angsd/get_maf/{population}.{chr}.depthSample`: depth histogram for
    each sample at SNP positions (see
    <https://www.popgen.dk/angsd/index.php/Allele_Counts> for details).
- Hard-coded arguments (see
  <https://www.popgen.dk/angsd/index.php/Allele_Frequencies> for
  details):
  - `-doGlf 2`: export beagle-formatted likelihood files.
  - `-doMaf 1`: return minor allele frequency at each site
  - `doSaf 1`: calculate the site allele frequency likelihoods for each
    population at each chromosome
  - `-sites {input.site_list}`: use the SNP list generated
    `snp_calling_global` rule to filter sites for analyses
  - `-doMajorMinor 3`: the major and minor allele are pre-determined
    based on the SNP list supplied to `-site`.
  - `-doCounts1`: return count of alleles at each site
  - `-dumpCounts1`: print the summed depth of all individuals at each
    site in the `.pos.gz` file
  - `-doIBS 1`: randomly print a single base from each individual at
    each position
  - `-makematrix 1`: print out a matrix of average distance between
    pairs of individuals
  - `-doCov 1`: print out the covariance matrix which can be used for
    PCA
  - `remove_bads 1`: remove duplicated or failed reads.
  - `-only_proper_pairs 1`: only include paired end data that both reads
    of a pair are mapped correctly
  - `-r {wildcards.chr}`: constrain the analysis on each chromosome
    separately. This ensures that all chromosomes are analyzed in
    parallel.
- Customizable arguments (see
  <https://www.popgen.dk/angsd/index.php/Allele_Frequencies> for
  details):
  - `-anc` : path to the reference genome. Users can specify their
    reference genome in the config file.
  - `-P` : number of threads this rule uses. This is defaulted to be 8
    but can be modified in the config file.
  - `-GL`: model chosen for genotype likelihood estimation. 1: SAMtools
    model, 2: GATK models, 3: SOAPsnp model, 4: SYK model.
  - `-minInd`: minimum number of individuals that have to have data for
    a site to be included. This is defaulted to be 1 but can be modified
    in the config file.
  - `-setMinDepthInd`: minimum read depth an individual must have to be
    included in the count of individuals for `-minInd`. This is
    defaulted to be 1 and can be modified in the config file.

## Fst estimation (`get_fst.smk`)

This Snakefile estimates Fst between pairs of populations and visualizes
the result.

#### get_fst

This rule estimates Fst between each population pair using the realSFS
module in ANGSD. It first uses the site allele frequency likelihoods
(saf) files of both populations to construct a chromosome-wide
two-dimensional site frequency spectrum (2dSFS). It then takes this 2D
SFS as a prior for the calculation of posterior estimates of Fst per
SNP. Lastly, it generates an average Fst estimate per chromosome. Note
that users need to specify whether their reference genome can represent
the ancestral state in the config file (e.g. if it comes from an
outgroup). If so, an unfolded SFS will be estimated. Otherwise, a folded
SFS will be estimated.

- Key output files
  - `angsd/get_fst/{population1}.{population2}.{chr}.2dSFS`: 2
    dimensional site frequency spectrum that will be used as a prior to
    estimate Fst (see
    <https://www.popgen.dk/angsd/index.php/2d_SFS_Estimation> for
    details).
  - `angsd/get_fst/{population1}.{population2}.{chr}.alpha_beta.txt`:
    variance partitioning at each SNP location used to calculate Fst.
    Weighted Fst in a window can be obtained by dividing the sum of the
    third column (`A`) by the sum of the fourth column (`B`) in that
    window (see
    <https://www.popgen.dk/angsd/index.php/Fst#realSFS_fst_print> for
    details).
  - `angsd/get_fst/{population1}.{population2}.{chr}.average_fst.txt`:
    estimated unweighted (the first value in this file) and weighted
    (the second value in this file) average Fst for each chromosome. The
    weighted average is preferred in most cases.
- Customizable arguments
  - `-P`: number of threads this rule uses. This is defaulted to be 8
    but can be modified in the config file.
  - `-fold`: folded vs. unfolded 2dSFS. An unfolded 2dSFS is estimated
    by default, but if users specify in the config file that the
    reference genome cannot represent the ancestral state (through the
    `ref_type` flag), a folded 2dSFS will be estimated instead.

#### plot_fst

This rule outputs three types of Manhattan plots for each pair of
populations: 1) per-SNP Fst, 2) Fst in windows with a fixed length, and
3) Fst in windows with a fixed number of SNPs.

- Key output files
  - `figures/fst/{population1}.{population2}.{snp_window_size}snp_window.png`:
    plot of fixed-SNP sliding window Fst
  - `figures/fst/{population1}.{population2}.{bp_window_size}bp_window.png`:
    plot of fixed-base pair sliding window Fst
  - `figures/fst/{population1}.{population2}.png`: plot of original Fst
    with no sliding window
- Customizable arguments
  - `snp_window_size`: the size of a fixed-SNP sliding window. The
    default value is 100 but it can be modified in the config file.
  - `bp_window_size`: the size of a fixed-base pair sliding window (in
    bp). The default value is 10000bp but it can be modified in the
    config file.
