rule get_theta:
    input:
        ref=REFERENCE,
        fai=REFERENCE+".fai",
        bamlist="{basedir}/docs/{population}_" + BAMLIST,
        site_list="{basedir}/angsd/get_depth_global/{chr}.site_list",
    output:
        expand("{{basedir}}/angsd/get_theta/{{population}}.{{chr}}.{ext}",
               ext=[ "arg", "sfs", "thetas.idx", "thetas.gz", "thetas.tsv.gz", "saf.gz", "saf.idx", "saf.pos.gz", "average_thetas.pestPG" ]),
        touch("{basedir}/angsd/get_theta/{population}.{chr}.done"),
    threads:
        config["get_theta"]["threads"]
    params:
        outdir="{basedir}/angsd/get_theta",
        outbase="{basedir}/angsd/get_theta/{population}.{chr}",
        minq=config["get_theta"]["minq"],
        minmapq=config["get_theta"]["minmapq"],
        window_size=config["get_theta"]["window_size"],
        step_size=config["get_theta"]["step_size"],
        dosaf_extra=config["get_theta"]["dosaf_extra"],
        realsfs_extra=config["get_theta"]["realsfs_extra"],
    log: "{basedir}/angsd/get_theta/{population}.{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        mkdir -p {params.outdir}
        ## Get saf file
        angsd \
        -b {input.bamlist} \
        -sites {input.site_list} -r {wildcards.chr} \
        -anc {input.ref} \
        -out {params.outbase} \
        -doSaf 1 \
        -doCounts 1 \
        -GL 1 \
        -P {threads} \
        -minQ {params.minq} -minMapQ {params.minmapq} \
        -remove_bads 1 -only_proper_pairs 1\
        {params.dosaf_extra} 2> {log}
        
        ## Get SFS from saf
        realSFS \
        {params.outbase}.saf.idx \
        -P {threads} \
        {params.realsfs_extra} \
        > {params.outbase}.sfs 2>> {log}
        
        ## Estimate theta
        realSFS saf2theta \
        {params.outbase}.saf.idx \
        -outname {params.outbase} \
        -sfs {params.outbase}.sfs \
        -anc {input.ref} \
        -P {threads} 2>> {log}
        
        ## Print per-SNP theta
        thetaStat print \
        {params.outbase}.thetas.idx | \
        gzip \
        > {params.outbase}.thetas.tsv.gz 2>> {log}
        
        ## Calculate per-chromosome average theta
        thetaStat do_stat \
        {params.outbase}.thetas.idx \
        -outnames {params.outbase}.average_thetas 2>> {log}
        
        ## Calculate fixed window theta
        thetaStat do_stat \
        {params.outbase}.thetas.idx \
        -win {params.window_size} -step {params.step_size} \
        -outnames {params.outbase}.{params.window_size}window_{params.step_size}step.thetas 2>> {log}
        '''

rule plot_theta_by_window:
    input:
        chr_table = CHR_TABLE,
        done = expand("{{basedir}}/angsd/get_theta/{{population}}.{chr}.done", chr = CHRS),
    output:
        plot = "{basedir}/figures/get_theta/{population}.theta_by_window.png",
        done = touch("{basedir}/angsd/get_theta/{population}.plot_theta_by_window.done"),
    params:
        indir = "{basedir}/angsd/get_theta",
        outdir = "{basedir}/figures/get_theta",
        window_size=config["get_theta"]["window_size"],
        step_size=config["get_theta"]["step_size"],
        rscript = config["global"]["scriptdir"] + "/plot_theta_by_window.R",
    threads: 4
    log: "{basedir}/angsd/get_theta/{population}.plot_theta_by_window.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {params.indir} {output.plot} {input.chr_table} {params.window_size} {params.step_size} {wildcards.population} &> {log}
        '''