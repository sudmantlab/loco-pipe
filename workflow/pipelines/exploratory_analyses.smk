import numpy as np
import pandas as pd
import os 

configfile: "config/exploratory_analyses.yaml"

#################################################
##  Things read directly from the config file

BASEDIR  = config["global"]["basedir"]
REFERENCE = config["global"]["reference"]
CHR_TABLE = config["global"]["chr_table"]
META_DATA = config["global"]["meta_data"]
BAMLIST = config["global"]["bamlist"]
POP_L1 = config["global"]["pop_level"] if config["to_include"]["run_pcangsd_local"]|config["to_include"]["get_maf"]| config["to_include"]["get_theta"]|config["to_include"]["get_Fst"] else []
#################################################

#################################################
## Varibles generated from the config variables

CHR_TB = pd.read_csv(CHR_TABLE, sep="\t",header=None, index_col=None)
CHRS = CHR_TB.iloc[:,0]
CHR_LIST_PATHWAY = BASEDIR + "/docs/chr_list.txt"
CHR_LIST = CHRS.to_csv(CHR_LIST_PATHWAY, header=False,index=False,sep="\t")
MD_TABLE = pd.read_csv(META_DATA, sep="\t")
L1_COL = MD_TABLE[POP_L1]
ALL_BAMS = MD_TABLE["bam"].unique()
ALL_SAMPLES = MD_TABLE["sample_name"].unique()
BAMERORY = MD_TABLE.set_index('sample_name')['bam'].to_dict()
df = pd.read_csv(META_DATA, sep="\t", usecols=[POP_L1]).sort_values(by=[POP_L1]).drop_duplicates()
GROUPS = df[POP_L1].unique()
#################################################

#################################################
## Build pathways for output files

POP_L1_PCANGSD_PATH = L1_COL.apply(lambda val: config["global"]["basedir"] + "/pcangsd/local/" + val).unique()
POP_L1_MAF_PATH  = L1_COL.apply(lambda val: config["global"]["basedir"] + "/angsd/get_maf/" + val).unique()
POP_L1_THETA_PATH =L1_COL.apply(lambda val: config["global"]["basedir"] + "/angsd/get_theta/" + val).unique()

POP_L1_FST_PLOT_PATH = []
POP_L1_FST_PATH = []
for i in range(len(GROUPS)):
  for j in range(i+1, len(GROUPS)):
    POP_L1_FST_PLOT_PATH.append(config["global"]["basedir"] + '/figures/fst/' + GROUPS[i] + '.' + GROUPS[j])
    POP_L1_FST_PATH.append(config["global"]["basedir"] + '/angsd/get_fst/' + GROUPS[i] + '.' + GROUPS[j])
#################################################

wildcard_constraints:
    chr = '|'.join([x for x in CHRS]),
    id =  '|'.join([x for x in ALL_SAMPLES]),
    
ruleorder: 
    snp_calling_global > subset_snp_list_global > combine_snp_list_global > combine_subsetted_beagle_global
ruleorder: 
    get_site_list_global > combine_site_list_global

#################################################
rule all:
    input:
        expand("{basedir}/docs/{bamlist}", basedir = BASEDIR,bamlist = BAMLIST),
        expand("{basedir}/docs/{population}_{bamlist}", basedir = BASEDIR, population = GROUPS, bamlist=BAMLIST),
        expand("{bam}.bai", bam=ALL_BAMS),
        expand("{basedir}/angsd/get_depth_global/{chr}.done", basedir = BASEDIR, chr = CHRS),
        expand("{basedir}/angsd/get_depth_global/depth_filter.tsv", basedir=BASEDIR),
        expand("{basedir}/angsd/get_depth_global/depth_filter.png", basedir=BASEDIR),
        expand("{basedir}/angsd/snp_calling_global/{chr}.done", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/heterozygosity/{id}.done", basedir=BASEDIR, id=ALL_SAMPLES),
        expand("{basedir}/angsd/snp_calling_global/combined.subsetted.beagle.gz", basedir=BASEDIR),
        
        expand("{basedir}/ohana/global/combined.subsetted.k{k}.done", 
               basedir=BASEDIR, 
               k = list(range(config["run_ohana_global"]["min_k"],
               config["run_ohana_global"]["max_k"] + 1))) if config["to_include"]["run_ohana_global"] else [],
               

        expand("{basedir}/ohana/local/{population}.combined.subsetted.k{k}.done", 
               basedir=BASEDIR, population = GROUPS,
               k = list(range(config["run_ohana_local"]["min_k"],
               config["run_ohana_local"]["max_k"] + 1))) if config["to_include"]["run_ohana_local"] else [],
               
        expand("{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list.done", basedir=BASEDIR),
        expand("{basedir}/pcangsd/global/combined.subsetted.done", basedir = BASEDIR) if config["to_include"]["run_pcangsd_global"] else [],
        expand("{path}.combined.subsetted.done", path = POP_L1_PCANGSD_PATH) if config["to_include"]["run_pcangsd_local"] else [],
        expand("{path}.{chr}.done", path=POP_L1_MAF_PATH, chr=CHRS) if config["to_include"]["get_maf"] else [],
        expand("{path}.{chr}.done", path=POP_L1_THETA_PATH, chr=CHRS) if config["to_include"]["get_theta"] else [],
        expand("{path}.{chr}.done", path=POP_L1_FST_PATH, chr=CHRS) if config["to_include"]["get_Fst"] else [],
        expand("{path}.done", path=POP_L1_FST_PLOT_PATH) if config["to_include"]["plot_Fst"] else [],




include: "../rules/get_bamlist_global.smk"
include: "../rules/get_bamlist_local.smk"
include: "../rules/index_bam.smk"
include: "../rules/get_depth_global.smk"
include: "../rules/get_depth_filter_global.smk"
include: "../rules/get_site_list_global.smk"
include: "../rules/get_heterozygosity.smk"
include: "../rules/snp_calling_global.smk"
include: "../rules/subset_snp_list.smk"
include: "../rules/combine_snp_list.smk"
include: "../rules/run_ohana.smk"
include: "../rules/get_fst.smk"
include: "../rules/run_pcangsd.smk"
include: "../rules/get_maf.smk"
include: "../rules/get_theta.smk"
