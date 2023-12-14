# This rule takes in population-specific bamlists and then estimate allele frequencies at each SNP location without subsetting. It will also output population-specific genotype likelihood in Beagle format and site allele frequency liklihoods.
rule get_maf:
    input:
        ref= REFERENCE,
        fai= REFERENCE+".fai",
        bamlist="{basedir}/docs/{population}_" + BAMLIST,
        site_list="{basedir}/angsd/snp_calling_global/{chr}.snp_list"
    output:
        expand("{{basedir}}/angsd/get_maf/{{population}}.{{chr}}.{ext}",
               ext=[ "arg", "beagle.gz", "covMat", "ibs.gz", "ibsMat", "mafs.gz", "pos.gz", "saf.gz", "saf.idx", "saf.pos.gz" ]),
        touch("{basedir}/angsd/get_maf/{population}.{chr}.done")
    threads:
        config["get_depth_global"]["threads"]
    params:
        outdir="{basedir}/angsd/get_maf",
        gl_model=config["global"]["gl_model"],
        minq=config["get_depth_global"]["minq"],
        minmapq=config["get_depth_global"]["minmapq"],
        extra=config["get_maf"]["extra"]
    log: "{basedir}/angsd/get_maf/{population}.{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        mkdir -p {params.outdir}
        angsd \
            -bam {input.bamlist} -anc {input.ref} \
            -P {threads} \
            -out {params.outdir}/{wildcards.population}.{wildcards.chr} \
            -sites {input.site_list} -r {wildcards.chr} \
            -GL {params.gl_model} -doGlf 2 -doMaf 1 -doSaf 1 -doMajorMinor 3 \
            -doCounts 1 -dumpCounts 1 \
            -doIBS 1 -makematrix 1 -doCov 1 \
            -minQ {params.minq} -minMapQ {params.minmapq} \
            -remove_bads 1 -only_proper_pairs 1 \
            {params.extra} 2> {log}
        '''