# This function returns all the index files for the input bam files.
# It is used in the "input" of the following rule.
def get_input_bai(lists):
    bai = expand("{list}.bai", list = lists)
    return bai

# This rule involves two separate steps from ANGSD. First, it accomplishes SNP Calling with the bam files,
# and then uses the "sites" method to output the snp indices for each chromosome (i.e. snp_list.idx).
# Detailed description on SNP Calling in ANGSD could be found on this webpage: http://www.popgen.dk/angsd/index.php/SNP_calling.
# Detailed description on "sites" method in ANGSD could be found on this webpage: http://www.popgen.dk/angsd/index.php/Sites.
rule snp_calling_global:
    input:
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        bamlist = "{basedir}/docs/" + BAMLIST,
        bai = get_input_bai(ALL_BAMS),
        depth_filter = "{basedir}/angsd/get_depth_global/depth_filter.tsv"
    output:
        expand("{{basedir}}/angsd/snp_calling_global/{{chr}}.{ext}",
               ext = [ "arg", "beagle.gz", "covMat", "depthGlobal", "depthSample", "ibs.gz", "ibsMat", "mafs.gz", "pos.gz", "snp_list", "snp_list.idx", "snp_list.bin" ]),
        touch("{basedir}/angsd/snp_calling_global/{chr}.done")
    threads:
        config["get_depth_global"]["threads"]
    params:
        outdir = "{basedir}/angsd/snp_calling_global",
        gl_model = config["global"]["gl_model"],
        minq = config["get_depth_global"]["minq"],
        minmapq = config["get_depth_global"]["minmapq"],
        minind = config["snp_calling_global"]["minind"],
        minmaf = config["snp_calling_global"]["minmaf"],
        pval = config["snp_calling_global"]["pval"],
        extra = config["snp_calling_global"]["extra"]
    log: "{basedir}/angsd/snp_calling_global/{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        MINDP=`cat {input.depth_filter} | tail -n 1 | cut -f 1`
        MAXDP=`cat {input.depth_filter} | tail -n 1 | cut -f 2`
        NIND=`cat {input.bamlist} | wc -l`
        MININD_PROPORTION={params.minind}
        MININD=`awk -v x="$NIND" -v y="$MININD_PROPORTION" 'BEGIN {{ z = int(x * y); print z }}'`
        mkdir -p {params.outdir}
        angsd \
            -bam {input.bamlist} -ref {input.ref} \
            -P {threads} \
            -out {params.outdir}/{wildcards.chr} -r {wildcards.chr} \
            -GL {params.gl_model} -doGlf 2 -doMaf 1 -doMajorMinor 1 \
            -doDepth 1 -doCounts 1 -maxDepth 100000 -dumpCounts 1 \
            -doIBS 1 -makematrix 1 -doCov 1 \
            -setMinDepth $MINDP -setMaxDepth $MAXDP -minInd $MININD \
            -SNP_pval {params.pval} -minMaf {params.minmaf} \
            -minQ {params.minq} -minMapQ {params.minmapq} \
            -remove_bads 1 -only_proper_pairs 1 \
            {params.extra} 2> {log}
        gunzip -c {params.outdir}/{wildcards.chr}.mafs.gz | cut -f 1,2,3,4 | tail -n +2 \
        > {params.outdir}/{wildcards.chr}.snp_list
        angsd sites index {params.outdir}/{wildcards.chr}.snp_list
        '''
