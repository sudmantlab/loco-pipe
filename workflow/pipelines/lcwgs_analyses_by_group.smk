import numpy as np
import pandas as pd
from functools import partial
import os 

configfile: "config/lcwgs_analyses_by_group.json"

META_DATA=config["global"]["meta_data_filtered"]
GROUPS=pd.read_csv(META_DATA, sep="\t")["species_group"].unique()
COMMON_NAMES=pd.read_csv(META_DATA, sep="\t")["common_name"].unique()
BASEDIR=expand(config["global"]["basedir"] + "/by_group/{group}", group=GROUPS)
REFERENCE=config["global"]["reference"]
CHR_LIST=config["global"]["chr_list"]
CHRS=pd.read_csv(CHR_LIST, header=None)[0].unique()
BAMLIST=config["get_depth_global"]["bamlist"]
SPECIES_LIST=config["global"]["species_list"]

#ALL_BAMS=pd.read_csv(BAMLIST, header=None)[0].unique()
#ALL_SAMPLES=[os.path.basename(bam).split('.')[0] for bam in ALL_BAMS]
#BAMDIR=config["global"]["bamdir"]

## get output path for maf estimation
SPECIES_LEVEL_MAF_PATH = pd.read_csv(META_DATA, sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/angsd/get_maf/' + x.common_name).path.unique()
PCA_CLUSTER_LEVEL_MAF_PATH = pd.read_csv(config["pca_cluster"]["pca_cluster_summary"], sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/angsd/get_maf/' + x.pca_cluster).path.unique()

## get output path for theta estimation
SPECIES_LEVEL_THETA_PATH = pd.read_csv(META_DATA, sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/angsd/get_theta/' + x.common_name).path.unique()

## get output path for fst estimation
df = pd.read_csv(META_DATA, sep="\t", usecols=['species_group', 'common_name']).sort_values(by=['species_group', 'common_name']).drop_duplicates()
groups = df.groupby('species_group')
SPECIES_LEVEL_FST_PLOT_PATH = []
SPECIES_LEVEL_FST_PATH = []
for name, group in groups:
    common_names = group['common_name'].tolist()
    for i in range(len(common_names)):
        for j in range(i+1, len(common_names)):
                SPECIES_LEVEL_FST_PLOT_PATH.append(config["global"]["basedir"] + '/figures/fst/' + common_names[i] + '.' + common_names[j])
                SPECIES_LEVEL_FST_PATH.append(config["global"]["basedir"] + '/angsd/get_fst/' + common_names[i] + '.' + common_names[j])
PCA_CLUSTER_LEVEL_FST_PLOT_PATH = pd.read_csv(config["pca_cluster"]["pca_cluster_pairs"], sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/figures/fst/' + x.pca_cluster_1 + '.' + x.pca_cluster_2).path.unique()

## get output path for species-level ohana analysis
SPECIES_LEVEL_OHANA_PATH = pd.read_csv(META_DATA, sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/ohana/local/' + x.common_name).path.unique()

## get output path for species-level ohana analysis
SPECIES_LEVEL_PCANGSD_PATH = pd.read_csv(META_DATA, sep="\t").assign(path=lambda x: config["global"]["basedir"] + '/by_group/' + x.species_group + '/pcangsd/local/' + x.common_name).path.unique()

wildcard_constraints:
    chr="\w+"
ruleorder: 
    snp_calling_global > subset_snp_list_global > combine_snp_list_global
ruleorder: 
    get_site_list_global > combine_site_list_global
ruleorder: 
    subset_beagle_global > combine_subsetted_beagle_global

rule all:
    input:
        expand("{basedir}/angsd/get_depth_global/{chr}.done", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/get_depth_global/depth_filter.tsv", basedir=BASEDIR),
        expand("{basedir}/angsd/get_depth_global/depth_filter.png", basedir=BASEDIR),
        expand("{basedir}/angsd/snp_calling_global/{chr}.done", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/snp_calling_global/{chr}.done", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/get_depth_global/{chr}.site_list.idx", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/snp_calling_global/combined.subsetted.beagle.gz", basedir=BASEDIR),
        expand("{basedir}/ohana/global/combined.subsetted.k{k}.q.matrix", 
               basedir=BASEDIR, 
               k=list(range(config["run_ohana_global"]["min_k"],
               config["run_ohana_global"]["max_k"] + 1))),
        expand("{basedir}/pcangsd/global/combined.subsetted.done", basedir=BASEDIR),
        expand("{path}.{chr}.done", path=SPECIES_LEVEL_MAF_PATH, chr=CHRS),
        expand("{path}.{chr}.done", path=SPECIES_LEVEL_FST_PATH, chr=CHRS),
        expand("{path}.{chr}.done", path=SPECIES_LEVEL_THETA_PATH, chr=CHRS),
        expand("{path}.done", path=SPECIES_LEVEL_FST_PLOT_PATH),
        expand("{path}.combined.subsetted.k{k}.q.matrix", 
               path=SPECIES_LEVEL_OHANA_PATH, 
               k=list(range(config["run_ohana_local"]["min_k"],
               config["run_ohana_local"]["max_k"] + 1))),
        expand("{path}.combined.subsetted.done", path=SPECIES_LEVEL_PCANGSD_PATH),
        expand("{path}.{chr}.done", path=PCA_CLUSTER_LEVEL_MAF_PATH, chr=CHRS),
        expand("{path}.done", path=PCA_CLUSTER_LEVEL_FST_PLOT_PATH),

include: "../rules/get_depth_global.smk"
include: "../rules/get_depth_filter_global.smk"
include: "../rules/get_site_list_global.smk"
include: "../rules/snp_calling_global.smk"
include: "../rules/subset_snp_list.smk"
include: "../rules/combine_snp_list.smk"
include: "../rules/run_ohana.smk"
include: "../rules/run_pcangsd.smk"
include: "../rules/get_maf.smk"
include: "../rules/get_fst.smk"
include: "../rules/get_theta.smk"

