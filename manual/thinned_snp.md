Analyses using the thinned SNP list (global or local)
================

- [Subset SNP list
  (`subset_snp_list.smk`)](#subset-snp-list-subset_snp_listsmk)
  - [subset_snp_list_global](#subset_snp_list_global)
  - [subset_beagle_global](#subset_beagle_global)
  - [subset_beagle_local](#subset_beagle_local)
- [Combine SNP list
  (`combine_snp_list.smk`)](#combine-snp-list-combine_snp_listsmk)
  - [combine_snp_list_global](#combine_snp_list_global)
  - [combine_subsetted_beagle_global](#combine_subsetted_beagle_global)
  - [combine_subsetted_beagle_local](#combine_subsetted_beagle_local)
- [PCA (`run_pcangsd.smk`)](#pca-run_pcangsdsmk)
  - [run_pcangsd_pca_global](#run_pcangsd_pca_global)
  - [plot_pcangsd_pca_global](#plot_pcangsd_pca_global)
  - [run_pcangsd_pca_local](#run_pcangsd_pca_local)
  - [plot_pcangsd_pca_local](#plot_pcangsd_pca_local)
- [Admixture analysis
  (`run_ohana.smk`)](#admixture-analysis-run_ohanasmk)
  - [get_ohana_input_global](#get_ohana_input_global)
  - [run_ohana_global](#run_ohana_global)
  - [plot_ohana_admixture_global](#plot_ohana_admixture_global)
  - [get_ohana_input_local](#get_ohana_input_local)
  - [run_ohana_local](#run_ohana_local)
  - [plot_ohana_admixture_local](#plot_ohana_admixture_local)

## Subset SNP list (`subset_snp_list.smk`)

This Snakefile prepares the thinned SNP lists and thinned genotype
likelihood files in beagle format needed for the analyses that require
unlinked markers.

#### subset_snp_list_global

This rule uses the global SNP lists outputted from the
`snp_calling_global` rule and generates thinned SNP lists by choosing
one in every “n” SNPs. This rule runs on each chromosome separately.

Customizable arguments n: thinning the SNP list by choosing one in every
“n” SNPs. It is defaulted to be 3 but can be modified in the config
file.

#### subset_beagle_global

This rule uses the global beagle formatted genotype likelihood files
outputted from the `snp_calling_global` rule and generates thinned
beagle files by choosing one in every “n” SNPs. This rule runs on each
chromosome separately. Customizable arguments n: thinning the SNP list
by choosing one in every “n” SNPs. It is defaulted to be 3 but can be
modified in the config file.

#### subset_beagle_local

This rule uses the population-specific beagle formatted genotype
likelihood files outputted from the `get_maf` rule and generates thinned
beagle files by choosing one in every “n” SNPs. This rule runs on each
chromosome separately. Customizable arguments n: thinning the SNP list
by choosing one in every “n” SNPs. It is defaulted to be 3 but can be
modified in the config file.

## Combine SNP list (`combine_snp_list.smk`)

This Snakefile combines SNP lists and beagle formatted genotype
likelihood files of individual chromosomes into a single file for
downstream analysis.

#### combine_snp_list_global

This rule combines the global full and thinned SNP lists of individual
chromosomes together into a single file. The resulting SNP lists are not
currently used for any downstream analyses (since the beagle files can
be directly passed to those analyses instead), but users may find them
helpful if they need to perform additional analyses at these SNP
locations using ANGSD or other software.

#### combine_subsetted_beagle_global

This rule combines the global thinned beagle formatted genotype
likelihoods of individual chromosomes into a single file for downstream
analysis.

#### combine_subsetted_beagle_local

This rule combines the population-specific thinned beagle formatted
genotype likelihoods of individual chromosomes into a single file for
downstream analysis.

## PCA (`run_pcangsd.smk`)

This Snakefile uses PCAngsd to conduct PCA. Please see
<https://github.com/Rosemeis/pcangsd> for details.

#### run_pcangsd_pca_global

This rule uses PCAngsd to conduct PCA with all samples combined using
the global thinned beagle file as input.

- Key output files
  - `pcangsd/global/combined.subsetted.cov`: a covariance matrix of all
    samples. Performing eigendecomposition on this matrix will yield
    each sample’s projection on principal component axes.
- Hard-coded arguments
  - `-snp_weights`: output the SNP weights of the significant K
    eigenvectors
  - `-sites_save`: choose to save the kept sites after filtering
- Customizable arguments
  - `-minMaf`: minimum minor allele frequency threshold for a site to be
    considered as a SNP. It is defaulted to be 0.05 but can be modified
    in the config file.
  - `-threads`: number of threads this rule uses. It is defaulted to be
    8 but can be modified in the config file.

#### plot_pcangsd_pca_global

This rule visualizes the results from global PCA analysis. Samples are
colored by a variable of the user’s choosing.

- Key output files
  - `figures/pcangsd/global/combined.subsetted.png`: PCA plot of the
    first 8 PC axes. Samples are colored by a variable of the user’s
    choice.
- Customizable arguments
  - `color_by`: users can choose to color the points by a column in the
    metadata table. This can be population, but can also be another
    variable.

#### run_pcangsd_pca_local

This rule uses PCAngsd to conduct PCA within each population using the
population-specific thinned beagle files as input. Each population is
analyzed separately.

- Key output files
  - `pcangsd/local/{population}.combined.subsetted.cov`: a covariance
    matrix of all samples in a population. Performing eigendecomposition
    on this matrix will yield each sample’s projection on principal
    component axes.
- Hard-coded arguments
  - `-snp_weights`: output the SNP weights of the significant K
    eigenvectors
  - `-sites_save`: choose to save the kept sites after filtering
- Customizable arguments
  - `-minMaf`: minimum minor allele frequency threshold for a site to be
    considered as a SNP. It is defaulted to be 0.05 but can be modified
    in the config file.
  - `-threads`: number of threads this rule uses. It is defaulted to be
    8 but can be modified in the config file.

#### plot_pcangsd_pca_local

This rule visualizes the results from population-level PCA. Samples are
colored by a variable of the user’s choosing. ’ \* Key output files \*
`figures/pcangsd/local/{population}.combined.subsetted.png`: PCA plot of
the first 8 PC axes for each population. Samples are colored by their
population.

- Customizable arguments
  - `color_by`: users can choose to color the points by a column in the
    metadata table. This can be population, but all points will show up
    as the same color this way. It is sometimes useful to first run this
    rule and color points by populations, and then supply another column
    to the metadata table that reflects a more subtle structure (which
    may be revealed by the first iteration of PCA). and rerun this rul.

## Admixture analysis (`run_ohana.smk`)

This Snakefile runs Ohana to conduct admixture analysis. See
<https://github.com/jade-cheng/ohana> for details.

#### get_ohana_input_global

This rule converts the thinned beagle formatted genotype likelihood file
to an “.lgm” format file which is needed for Ohana.

- Key output files
  - `ohana/global/combined.subsetted.lgm`: genotype likelihoods in the
    format required by Ohana.

#### run_ohana_global

This rule performs admixture analysis with Ohana using the `qpas`
function. Analyses will be run across a range of k values (i.e. the
number of ancestral/source populations). Each k value will be run
separately.

- Key output files
  - `ohana/global/combined.subsetted.k{k}.q.matrix`: estimated admixture
    proportions of all samples for each value of k.
  - `ohana/global/combined.subsetted.k{k}.f.matrix`: estimated allele
    frequencies in source populations for each value of k.
- Hard-coded arguments:
  - `-mi 50`: maximum number of iterations is set to be 50 as
    recommended by Ohana’s author.
- Customizable arguments
  - `-k`: number of assumed ancestral/source populations. Its default
    range is 2-7 but can be modified in the config file.

#### plot_ohana_admixture_global

This rule generates an admixture plot. All k values as specified in the
config file will be included in the same plot.

- Key output files
  - `figures/ohana/global/combined.subsetted.png`: an admixture plot.
    Rows are arranged by different values of k, and columns are arranged
    based on the user’s specification.
- Customizable arguments
  - `group_by`: users can choose to group the samples by a column in the
    metadata table. This can be population, but can also be another
    variable.

#### get_ohana_input_local

This rule is the population-level counterpart of
`get_ohana_input_global`.

#### run_ohana_local

This rule is the population-level counterpart of `run_ohana_global`.

#### plot_ohana_admixture_local

This rule is the population-level counterpart of
`plot_ohana_admixture_global`.
