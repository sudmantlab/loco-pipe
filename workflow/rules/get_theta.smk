# This rule includes three steps.
# The first step is to input a population-specific bamlist and a chromosome-specific site list into ANGSD to calculate site allele frequency likelihoods for each population at each chromosome. Then we use realSFS to estimate the one-dimensional SFS on all sites (variable and invariable). Finally, the SFS is used as a prior to calculate theta. We offer three types of theta estimation: 1) per-SNP, 2) chromosome-average, and 3) sliding-window. 
rule get_theta:
    input:
        ref=REFERENCE,
        fai=REFERENCE+".fai",
        bamlist="{basedir}/docs/{population}_" + BAMLIST,
        site_list="{basedir}/angsd/get_depth_global/{chr}.site_list",
    output:
        expand("{{basedir}}/angsd/get_theta/{{population}}.{{chr}}.{ext}",
               ext=[ "arg", "sfs", "thetas.idx", "thetas.gz", "thetas.tsv.gz", "saf.gz", "saf.idx", "saf.pos.gz", "average_thetas.pestPG" ]),
        done = touch("{basedir}/angsd/get_theta/{population}.{chr}.done"),
    threads:
        config["get_theta"]["threads"]
    params:
        outdir="{basedir}/angsd/get_theta",
        outbase="{basedir}/angsd/get_theta/{population}.{chr}",
        gl_model=config["global"]["gl_model"],
        minq=config["get_theta"]["minq"],
        minind = config["get_theta"]["minind"],
        mindepthind = config["get_theta"]["setMinDepthInd"],
        minmapq=config["get_theta"]["minmapq"],
        ref_type = config["global"]["ref_type"],
        window_size=config["get_theta"]["window_size"],
        step_size=config["get_theta"]["step_size"],
        dosaf_extra=config["get_theta"]["dosaf_extra"],
        realsfs_extra=config["get_theta"]["realsfs_extra"],
    log: "{basedir}/angsd/get_theta/{population}.{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        MININD={params.minind}
        MINDEPTHIND={params.mindepthind}
        mkdir -p {params.outdir}
        ## Get saf file
        angsd \
        -b {input.bamlist} \
        -sites {input.site_list} -r {wildcards.chr} \
        -anc {input.ref} \
        -out {params.outbase} \
        -doSaf 1 \
        -doCounts 1 \
        -GL {params.gl_model} \
        -P {threads} \
        -minInd $MININD -setMinDepthInd $MINDEPTHIND \
        -minQ {params.minq} -minMapQ {params.minmapq} \
        -remove_bads 1 -only_proper_pairs 1\
        {params.dosaf_extra} 2> {log}
        
        ## Get SFS from saf
        realSFS \
        {params.outbase}.saf.idx \
        -P {threads} \
        -fold {params.ref_type} \
        {params.realsfs_extra} \
        > {params.outbase}.sfs 2>> {log}
      
        ## Estimate theta
        realSFS saf2theta \
        {params.outbase}.saf.idx \
        -outname {params.outbase} \
        -sfs {params.outbase}.sfs \
        -anc {input.ref} \
        -P {threads} \
        -fold {params.ref_type} \
        2>> {log}
        
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
# This rule generates a barplot for sfs distribution of each population.
rule plot_sfs_distribution:
    input:
        chr_list = "{basedir}/docs/chr_list.txt",
        done = expand("{{basedir}}/angsd/get_theta/{{population}}.{chr}.done", chr = CHRS),
    output:
        plot = "{basedir}/figures/theta/sfs_distribution/{population}.sfs.png",
        done = touch("{basedir}/figures/theta/sfs_distribution/{population}.sfs.done"),
    params:
        indir = "{basedir}/angsd/get_theta",
        outdir = "{basedir}/figures/theta/sfs_distribution",
        ref_type = config["global"]["ref_type"],
        rscript=config["global"]["scriptdir"] + "/plot_SFS.R",
    threads: 4
    log: "{basedir}/angsd/get_theta/sfs_distribution/{population}.plot_sfs_distribution.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript}  {params.indir} {params.outdir} {wildcards.population} {input.chr_list} {params.ref_type} &> {log}
        '''


       
# This rule generates plots for 1)estimates of pi, 2)Watterson’s theta, and 3)Tajima’s D in sliding windows for each population separately.
rule plot_theta_by_window:
    input:
        chr_table = CHR_TABLE_PATH,
        done = expand("{{basedir}}/angsd/get_theta/{{population}}.{chr}.done", chr = CHRS),
    output:
        plot = "{basedir}/figures/theta/{population}.theta_by_window.png",
        done = touch("{basedir}/figures/theta/{population}.theta_by_window.done"),
    params:
        indir = "{basedir}/angsd/get_theta",
        outdir = "{basedir}/figures/theta",
        window_size=config["get_theta"]["window_size"],
        step_size=config["get_theta"]["step_size"],
        fig_height = config["get_theta"]["fig_height"],
        fig_width = config["get_theta"]["fig_width"],
        rscript = config["global"]["scriptdir"] + "/plot_theta_by_window.R",
    threads: 4
    log: "{basedir}/angsd/get_theta/{population}.plot_theta_by_window.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {params.indir} {output.plot} {input.chr_table} {params.window_size} {params.step_size} {wildcards.population} {params.fig_height} {params.fig_width} &> {log}
        '''