# This function returns all the index files for the input bam files.
# It is used in the "input" of the following rule.
def get_input_bai(lists):
    bai = expand("{list}.bai", list = lists)
    return bai
    
# This rule runs ANGSD in a conda environment, counts the read depth of every site summed over all samples, 
# and outputs the result into a pos.gz file. It also outputs depthGlobal and depthSample,
# whose detailed description could be found on the ANGSD website:
# http://www.popgen.dk/angsd/index.php/Depth. 
rule get_depth_global:
    input:
        ref = REFERENCE,
        fai = REFERENCE + ".fai",
        bamlist = "{basedir}/docs/" + BAMLIST,
        bai = get_input_bai(ALL_BAMS),
    output:
        expand("{{basedir}}/angsd/get_depth_global/{{chr}}.{ext}",
               ext = [ "arg", "depthGlobal", "depthSample", "pos.gz" ]),
        touch("{basedir}/angsd/get_depth_global/{chr}.done")
    threads: 
        config["get_depth_global"]["threads"]
    params:
        outdir = "{basedir}/angsd/get_depth_global",
        minq = config["get_depth_global"]["minq"],
        minmapq = config["get_depth_global"]["minmapq"],
        extra = config["get_depth_global"]["extra"]
    log: "{basedir}/angsd/get_depth_global/{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        mkdir -p {params.outdir}
        angsd -bam {input.bamlist} -ref {input.ref} -P {threads} \
        -out {params.outdir}/{wildcards.chr} \
        -r {wildcards.chr} \
        -doDepth 1 -doCounts 1 -maxDepth 100000 -dumpCounts 1 \
        -minQ {params.minq} \
        -minMapQ {params.minmapq} \
        -remove_bads 1 -only_proper_pairs 1 \
        {params.extra} \
        2> {log}
        '''
