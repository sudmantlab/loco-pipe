# This rule uses PCAngsd software to conduct PCA analysis on a global level. 
# More details about PCAngsd could be found here: http://www.popgen.dk/software/index.php/PCAngsd.
rule run_pcangsd_pca_global:
    input: "{basedir}/angsd/snp_calling_global/{file}.beagle.gz",
    output: 
        cov = "{basedir}/pcangsd/global/{file}.cov",
        done = touch("{basedir}/pcangsd/global/{file}.done"),
    params:
        outdir = "{basedir}/pcangsd/global",
        minmaf=config["run_pcangsd_global"]["minmaf"],
    threads: config["run_pcangsd_global"]["threads"]
    log: "{basedir}/pcangsd/global/{file}.log"
    conda:
        "pcangsd_lcpipe"

    shell:
        '''
        conda env list > {log}
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp_weights \
        --sites_save \
        --minMaf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/{wildcards.file} \
        &>> {log}
        '''

# This rule uses PCAngsd software to conduct PCA analysis on a local level. 
# More details about PCAngsd could be found here: http://www.popgen.dk/software/index.php/PCAngsd.      
rule run_pcangsd_pca_local:
    input: "{basedir}/angsd/get_maf/{population}.{file}.beagle.gz",
    output: 
        cov = "{basedir}/pcangsd/local/{population}.{file}.cov",
        done = touch("{basedir}/pcangsd/local/{population}.{file}.done"),
    params:
        outdir = "{basedir}/pcangsd/local",
        minmaf=config["run_pcangsd_local"]["minmaf"],
    threads: config["run_pcangsd_local"]["threads"]
    log: "{basedir}/pcangsd/local/{population}.{file}.log"
    conda:
        "pcangsd_lcpipe"
    shell:
        '''
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp_weights \
        --sites_save \
        --minMaf {params.minmaf} \
        --threads {threads} \
        --out {params.outdir}/{wildcards.population}.{wildcards.file} \
        &> {log}
        '''
# This rule plots the PCA analysis conducted in the rule,run_pcangsd_pca_global, above. 
rule plot_pcangsd_pca_global:
    input:
        cov = "{basedir}/pcangsd/global/{file}.cov",
        done = "{basedir}/pcangsd/global/{file}.done",
    output:
        plot = "{basedir}/figures/pcangsd/global/{file}.png",
        done = touch("{basedir}/figures/pcangsd/global/{file}.done"),
    params:
        outdir = "{basedir}/figures/pcangsd/global",
        sample_table_path = "{basedir}/docs/" + config["global"]["sample_table"],
        color_by = config["run_pcangsd_global"]["color_by"],
        rscript = config["global"]["scriptdir"] + "/plot_pcangsd_pca.R",
    threads: 1
    log: "{basedir}/pcangsd/global/{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {input.cov} {output.plot} {params.sample_table_path} {params.color_by} &> {log}
        '''

# This rule plots the PCA analysis conducted in the rule,run_pcangsd_pca_local, above. 
rule plot_pcangsd_pca_local:
    input:
        cov = "{basedir}/pcangsd/local/{population}.{file}.cov",
        done = "{basedir}/pcangsd/local/{population}.{file}.done",
    output:
        plot = "{basedir}/figures/pcangsd/local/{population}.{file}.png",
        done = touch("{basedir}/figures/pcangsd/local/{population}.{file}.done"),
    params:
        outdir = "{basedir}/figures/pcangsd/local",
        sample_table_path = "{basedir}/docs/" + config["global"]["sample_table"],
        color_by = config["run_pcangsd_local"]["color_by"],
        pop_col = config["global"]["pop_level"],
        rscript = config["global"]["scriptdir"] + "/plot_pcangsd_pca.R",
    threads: 1
    log: "{basedir}/pcangsd/local/{population}.{file}.plot.log"
    conda:
        "../envs/r.yaml" 
    shell:
        '''
        mkdir -p {params.outdir}
        Rscript --vanilla {params.rscript} {input.cov} {output.plot} {params.sample_table_path} {params.color_by} {params.pop_col} {wildcards.population} &> {log}
        '''
