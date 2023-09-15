## Note that this assumes the saf is polarized with an ancestral genome. If this is not the case, use extra argument in the config file to output an folded SFS.

rule get_fst:
    input:
        ref = REFERENCE,
        fai = REFERENCE+".fai",
        saf1 = "{basedir}/angsd/get_maf/{population1}.{chr}.saf.idx",
        saf2 = "{basedir}/angsd/get_maf/{population2}.{chr}.saf.idx",
    output:
        expand("{{basedir}}/angsd/get_fst/{{population1}}.{{population2}}.{{chr}}.{ext}",
               ext = [ "2dSFS", "alpha_beta.txt", "alpha_beta.fst.idx", "alpha_beta.fst.gz", "average_fst.txt" ]),
        touch("{basedir}/angsd/get_fst/{population1}.{population2}.{chr}.done"),
    threads:
        config["get_fst"]["threads"]
    params:
        outdir = "{basedir}/angsd/get_fst",
        outbase = "{basedir}/angsd/get_fst/{population1}.{population2}.{chr}",
        extra = config["get_fst"]["extra"],
    log: "{basedir}/angsd/get_fst/{population1}.{population2}.{chr}.log"
    conda: "../envs/angsd.yaml"
    shell:
        '''
        mkdir -p {params.outdir}
        # Generate the 2dSFS to be used as a prior for Fst estimation (and individual plots)
        realSFS {input.saf1} {input.saf2} -P {threads} {params.extra} > {params.outbase}.2dSFS 2> {log}
        # Estimating Fst in angsd
        realSFS fst index  {input.saf1} {input.saf2} -sfs {params.outbase}.2dSFS -fstout {params.outbase}.alpha_beta
        realSFS fst print {params.outbase}.alpha_beta.fst.idx > {params.outbase}.alpha_beta.txt 2>> {log}
        # Estimating average Fst in angsd
        realSFS fst stats {params.outbase}.alpha_beta.fst.idx > {params.outbase}.average_fst.txt 2>> {log}
        '''

rule plot_fst:
    input: 
        fst = expand("{{basedir}}/angsd/get_fst/{{population1}}.{{population2}}.{chr}.done", chr = CHRS),
        chr_table = config["global"]["chr_table"],
    output: 
        png = "{basedir}/figures/fst/{population1}.{population2}.png",
        done = touch("{basedir}/figures/fst/{population1}.{population2}.done"),
    threads: 4
    params:
        rscript = config["global"]["scriptdir"] + "/plot_fst.R",
        snp_window_size = config["plot_fst"]["snp_window_size"],
        bp_window_size = config["plot_fst"]["bp_window_size"],
    log: "{basedir}/figures/fst/{population1}.{population2}.log"
    conda: "../envs/r.yaml"
    shell:
        '''
        # module load r  
        Rscript {params.rscript} {input.chr_table} {wildcards.basedir} {wildcards.population1} {wildcards.population2} {params.snp_window_size} {params.bp_window_size} &> {log}
        '''
