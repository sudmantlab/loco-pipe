# This rule operates on a two-round two-step mechanism. In the first round, it first combines every snp_list file of each chromosome into 
# one snp_list file, and then uses the "sites" method to output the snp indices for the combined snp_list file (i.e. combined_snp_list.idx).
# In the second round, it repeats the two steps on chromosome-specific subsetted.snp_lists and outputs indexed snp_list file
# (i.e. combined.subsetted.snp_list.idx).  
rule combine_snp_list_global:
    input: 
        snp_list = expand("{{basedir}}/angsd/snp_calling_global/{chr}.snp_list", chr = CHRS),
        subsetted_snp_list = expand("{{basedir}}/angsd/snp_calling_global/{chr}.subsetted.snp_list", chr = CHRS),
    output:
        snp_list = "{basedir}/angsd/snp_calling_global/combined.snp_list",
        snp_list_idx = "{basedir}/angsd/snp_calling_global/combined.snp_list.idx",
        snp_list_bin = "{basedir}/angsd/snp_calling_global/combined.snp_list.bin",
        subsetted_snp_list = "{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list",
        subsetted_snp_list_idx = "{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list.idx",
        subsetted_snp_list_bin = "{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list.bin",
        done = touch("{basedir}/angsd/snp_calling_global/combined.subsetted.snp_list.done"),
    threads: 1
    conda: "../envs/angsd.yaml"
    log: "{basedir}/angsd/snp_calling_global/combine_snp_list.log"
    shell:
        '''
        cat {input.snp_list} > {output.snp_list} 2> {log}
        angsd sites index {output.snp_list} 2>> {log}
        cat {input.subsetted_snp_list} > {output.subsetted_snp_list} 2>> {log}
        angsd sites index {output.subsetted_snp_list} 2>> {log}
        '''
        
# This rule extracts the header row of the first file, combines the non-header rows from the subsetted.beagle file for each
# chromosome into one file (i.e. subsetted.beagle), and then compresses it.
rule combine_subsetted_beagle_global:
    input: 
        subsetted_beagle = expand("{{basedir}}/angsd/snp_calling_global/{chr}.subsetted.beagle", chr = CHRS),
    output:
        subsetted_beagle = "{basedir}/angsd/snp_calling_global/combined.subsetted.beagle",
        subsetted_beagle_gz = "{basedir}/angsd/snp_calling_global/combined.subsetted.beagle.gz",
        done = touch("{basedir}/angsd/snp_calling_global/combined.subsetted.beagle.done"),
    threads: 1
    log: "{basedir}/angsd/snp_calling_global/combine_subsetted_beagle.log"
    shell:
        '''
        head -n 1 {input.subsetted_beagle[0]} > {output.subsetted_beagle} 2> {log}
        tail -q -n +2 {input.subsetted_beagle} >> {output.subsetted_beagle} 2>> {log}
        gzip -c {output.subsetted_beagle} > {output.subsetted_beagle}.gz 2>> {log}
        '''
        
# This rule applies the same logic as the one above (i.e. combine_subsetted_beagle_global), but instead of operating on the
# global genotype likelihood file in beagle format, it combines the population-specific beagle files.
rule combine_subsetted_beagle_local:
    input: 
        subsetted_beagle = expand("{{basedir}}/angsd/get_maf/{{population}}.{chr}.subsetted.beagle", chr=CHRS),
    output:
        subsetted_beagle = "{basedir}/angsd/get_maf/{population}.combined.subsetted.beagle",
        subsetted_beagle_site = "{basedir}/angsd/get_maf/{population}.combined.subsetted.beagle.sites",
        subsetted_beagle_gz = "{basedir}/angsd/get_maf/{population}.combined.subsetted.beagle.gz",
        done = touch("{basedir}/angsd/get_maf/{population}.combined.subsetted.beagle.done"),
    threads: 1
    log: "{basedir}/angsd/get_maf/{population}.combine_subsetted_beagle.log"
    shell:
        '''
        head -n 1 {input.subsetted_beagle[0]} > {output.subsetted_beagle} 2> {log}
        tail -q -n +2 {input.subsetted_beagle} >> {output.subsetted_beagle} 2>> {log}
        cut -f 1 {output.subsetted_beagle} > {output.subsetted_beagle_site} 2>> {log}
        gzip -c {output.subsetted_beagle} > {output.subsetted_beagle_gz} 2>> {log}
        '''

