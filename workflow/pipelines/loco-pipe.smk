import numpy as np
import pandas as pd
import os 

configfile: "config/loco-pipe.yaml"

#################################################
##  Things read directly from the config file

BASEDIR  = config["global"]["basedir"]
REFERENCE = config["global"]["reference"]
CHR_TABLE = BASEDIR + "/docs/" + config["global"]["chr_table"]
METADATA = BASEDIR + "/docs/" + config["global"]["metadata"]

## Note: we need to define a flag for all population-level analyses, and this flag needs to be consistently applied
POP_L1_COLNAME = config["global"]["pop_level"] if config["to_include"]["run_pcangsd_local"]|config["to_include"]["get_maf"]| config["to_include"]["get_theta"]|config["to_include"]["get_Fst"] else []

#################################################
## Varibles generated from the config variables

CHR_TB = pd.read_csv(CHR_TABLE, sep="\t",header=None, index_col=None)
CHRS = CHR_TB.iloc[:,0]
#CHR_LIST = ','.join([x for x in CHRS])
#CHRS.to_csv(CHR_LIST, header=False,index=False,sep="\t") if os.path.exists(CHR_LIST) else []
MD_TABLE = pd.read_csv(METADATA, sep="\t")
ALL_BAMS = MD_TABLE["bam"].unique()
ALL_SAMPLES = MD_TABLE["sample_name"].unique()
BAM_DICT = MD_TABLE.set_index('sample_name')['bam'].to_dict()
POP_L1_COL = MD_TABLE[POP_L1_COLNAME]
POP_L1 = POP_L1_COL.unique()

#################################################
## Other variables

BAMLIST = "bamlist.txt"

#################################################
## Build pathways for output files

#POP_L1_PCANGSD_PATH = POP_L1_COL.apply(lambda val: config["global"]["basedir"] + "/pcangsd/local/" + val).unique()
#POP_L1_MAF_PATH  = POP_L1_COL.apply(lambda val: config["global"]["basedir"] + "/angsd/get_maf/" + val).unique()
#POP_L1_THETA_PATH = POP_L1_COL.apply(lambda val: config["global"]["basedir"] + "/angsd/get_theta/" + val).unique()

POP_L1_FST_PLOT_PATH = []
POP_L1_FST_PATH = []
for i in range(len(POP_L1)):
  for j in range(i+1, len(POP_L1)):
    POP_L1_FST_PLOT_PATH.append(config["global"]["basedir"] + '/figures/fst/' + POP_L1[i] + '.' + POP_L1[j])
    POP_L1_FST_PATH.append(config["global"]["basedir"] + '/angsd/get_fst/' + POP_L1[i] + '.' + POP_L1[j])
#################################################

wildcard_constraints:
    chr = '|'.join([x for x in CHRS]),
    id =  '|'.join([x for x in ALL_SAMPLES]),
    population =  '|'.join([x for x in POP_L1]),
    k = '[0-9]+'
    
ruleorder: 
    snp_calling_global > subset_snp_list_global > combine_snp_list_global > combine_subsetted_beagle_global
ruleorder: 
    get_site_list_global > combine_site_list_global

#################################################
rule all:
    input:
        ## bamlists
        expand("{basedir}/docs/{bamlist}", basedir = BASEDIR,bamlist = BAMLIST),
        expand("{basedir}/docs/{population}_{bamlist}", basedir = BASEDIR, population = POP_L1, bamlist=BAMLIST),
        
        ## bam indices
        expand("{bam}.bai", bam=ALL_BAMS),
        
        ## get depth and depth filter
        expand("{basedir}/angsd/get_depth_global/{chr}.done", basedir = BASEDIR, chr = CHRS),
        expand("{basedir}/angsd/get_depth_global/depth_filter.tsv", basedir=BASEDIR),
        expand("{basedir}/figures/depth/depth_filter.png", basedir=BASEDIR),
        
        ## SNP calling
        expand("{basedir}/angsd/snp_calling_global/{chr}.done", basedir=BASEDIR, chr=CHRS),
        expand("{basedir}/angsd/snp_calling_global/combined.subsetted.beagle.gz", basedir=BASEDIR),
        
        ## heterozygosity
        expand("{basedir}/angsd/heterozygosity/{id}.done", basedir=BASEDIR, id=ALL_SAMPLES) if config["to_include"]["get_heterozygosity"] else [],
        expand("{basedir}/angsd/heterozygosity/plot_heterozygosity.done", basedir=BASEDIR) if config["to_include"]["get_heterozygosity"] else [],
        
        ## Ohana global
        expand("{basedir}/ohana/global/combined.subsetted.k{k}.done", 
               basedir=BASEDIR, 
               k = list(range(config["run_ohana_global"]["min_k"],
               config["run_ohana_global"]["max_k"] + 1))) if config["to_include"]["run_ohana_global"] else [],
        expand("{basedir}/ohana/global/combined.subsetted.plot.done", basedir = BASEDIR) if config["to_include"]["run_ohana_global"] else [],

        ## Ohana local
        expand("{basedir}/ohana/local/{population}.combined.subsetted.k{k}.done", 
               basedir=BASEDIR, population = POP_L1,
               k = list(range(config["run_ohana_local"]["min_k"],
               config["run_ohana_local"]["max_k"] + 1))) if config["to_include"]["run_ohana_local"] else [],
        expand("{basedir}/ohana/local/{population}.combined.subsetted.plot.done", basedir = BASEDIR, population = POP_L1) if config["to_include"]["run_ohana_global"] else [],
        
        ### combine and subset SNP list
        expand("{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list.done", basedir=BASEDIR),
        
        ## PCAngsd global
        expand("{basedir}/pcangsd/global/combined.subsetted.done", basedir = BASEDIR) if config["to_include"]["run_pcangsd_global"] else [],
        expand("{basedir}/pcangsd/global/combined.subsetted.plot.done", basedir = BASEDIR) if config["to_include"]["run_pcangsd_global"] else [],
        
        ## PCAngsd local
        expand("{basedir}/pcangsd/local/{population}.combined.subsetted.done", basedir=BASEDIR, population = POP_L1) if config["to_include"]["run_pcangsd_local"] else [],
        expand("{basedir}/pcangsd/local/{population}.combined.subsetted.plot.done", basedir=BASEDIR, population = POP_L1) if config["to_include"]["run_pcangsd_local"] else [],
        
        ## maf and population-level genotype likelihood estimation
        expand("{basedir}/angsd/get_maf/{population}.{chr}.done",  basedir=BASEDIR, population = POP_L1, chr=CHRS) if config["to_include"]["get_maf"] else [],
        expand("{basedir}/angsd/get_maf/{population}.combined.subsetted.beagle.gz",  basedir=BASEDIR, population = POP_L1, chr=CHRS) if config["to_include"]["get_maf"] else [],
        
        ## theta and neutrality stats
        expand("{basedir}/angsd/get_theta/{population}.{chr}.done", basedir=BASEDIR, population = POP_L1, chr=CHRS) if config["to_include"]["get_theta"] else [],
        expand("{basedir}/angsd/get_theta/{population}.plot_theta_by_window.done", basedir=BASEDIR, population = POP_L1) if config["to_include"]["get_theta"] else [],

        ## Fst
        expand("{path}.{chr}.done", path=POP_L1_FST_PATH, chr=CHRS) if config["to_include"]["get_Fst"] else [],
        expand("{path}.done", path=POP_L1_FST_PLOT_PATH) if config["to_include"]["plot_Fst"] else [],
        
        ## local PCA with lostruct
        expand("{basedir}/lostruct/global/plot_lostruct_outlier_pca/plot_lostruct_outlier_pca.done", basedir = BASEDIR) if config["to_include"]["run_lostruct_global"] else [],




include: "../rules/pipeline_prep.smk"
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
include: "../rules/run_lostruct_global.smk"
