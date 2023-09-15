# This rule extracts the min and the max depth filters from the depth_filter.tsv first. With these filters, it subsets the pos.gz
# file generated from get_depth_global.smk for every chromosome. This generates a list of sites (both variable and invariable ones).
# Through the "-sites" flag in ANGSD, we can limit future analysis (e.g. genetic diveristy estimation) to this site list.
# A detailed explanation on "-sites" flag could be found on this webpage: http://www.popgen.dk/angsd/index.php/Sites.
rule get_site_list_global:
    input: 
        depth_filter = "{basedir}/angsd/get_depth_global/depth_filter.tsv",
        pos = "{basedir}/angsd/get_depth_global/{chr}.pos.gz",
    output:
        site_list = "{basedir}/angsd/get_depth_global/{chr}.site_list",
        site_list_idx = "{basedir}/angsd/get_depth_global/{chr}.site_list.idx"
    conda: "../envs/angsd.yaml"
    threads: 1
    shell:
        '''
        MINDP=`cat {input.depth_filter} | tail -n 1 | cut -f 3`
        MAXDP=`cat {input.depth_filter} | tail -n 1 | cut -f 4`
        zcat {input.pos} | awk -v min="$MINDP" -v max="$MAXDP" -F'\t' '$3>min && $3<max {{print $1"\t"$2}}' > {output.site_list}
        # Generate an index file for the site list.
        angsd sites index {output.site_list} 
        '''
        
# This rule collects the outputs from the previous rule (i.e. get_site_list_global) and combines filtered position
# of all chromosomes into one file. Then, similar to above, ANGSD is used to index the site list.
rule combine_site_list_global:
    input: 
        expand("{{basedir}}/angsd/get_depth_global/{chr}.site_list", chr = CHRS)
    output:
        site_list = "{basedir}/angsd/get_depth_global/combined.site_list",
        site_list_idx = "{basedir}/angsd/get_depth_global/combined.site_list.idx"
    threads: 1
    conda: "../envs/angsd.yaml"
    shell:
        '''
        cat {input} > {output.site_list}
        angsd sites index {output.site_list}
        '''
      
