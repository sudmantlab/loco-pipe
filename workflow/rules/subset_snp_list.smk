# This rule uses the snp_list file outputted from the snp_calling_global rule and generates a subset of
# SNPs for each chromosome by choosing one in every "n" SNPs. The resulting SNP list would be useful for
# analyses that require unlinked SNPs (e.g. admixture analysis).
rule subset_snp_list_global:
    input: 
        snp_list = "{basedir}/angsd/snp_calling_global/{chr}.snp_list",
    output: 
        subsetted_snp_list = "{basedir}/angsd/snp_calling_global/{chr}.subsetted.snp_list",
    params: n = config["subset_snp_list_global"]["n"]
    threads: 1
    shell:
        '''
        awk -v n={params.n} 'NR%n==0' {input.snp_list} > {output.subsetted_snp_list}
        '''
        
# This rule generates a subsetted beagle-formatted genotype likelihood file for each chromosome with a similar logic from 
# "subset_snp_list_global". Again, these files would be used for analyses that require unlinked SNPs. 
rule subset_beagle_global:
    input: 
        beagle = "{basedir}/angsd/snp_calling_global/{chr}.beagle.gz"
    output: 
        subsetted_beagle = "{basedir}/angsd/snp_calling_global/{chr}.subsetted.beagle",
    params: n = config["subset_snp_list_global"]["n"]
    threads: 1
    shell:
        '''
        zcat {input.beagle} | awk -v n={params.n} 'NR%n==1' > {output.subsetted_beagle}
        '''

# This rule generates a subset of beagle-formatted genotype likelihood file for each chromosome from one specific population. 
# The logic used to subset is the same as above.
rule subset_beagle_local:
    input: 
        beagle = "{basedir}/angsd/get_maf/{population}.{chr}.beagle.gz"
    output: 
        subsetted_beagle = "{basedir}/angsd/get_maf/{population}.{chr}.subsetted.beagle",
    params: n = config["subset_snp_list_local"]["n"]
    threads: 1
    shell:
        '''
        zcat {input.beagle} | awk -v n={params.n} 'NR%n==1' > {output.subsetted_beagle}
        '''
