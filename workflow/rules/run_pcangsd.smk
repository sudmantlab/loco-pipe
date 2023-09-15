rule run_pcangsd_global:
    input: "{basedir}/angsd/snp_calling_global/{file}.beagle.gz",
    output: 
        done = touch("{basedir}/pcangsd/global/{file}.done"),
    params:
        outdir = "{basedir}/pcangsd/global",
    threads: 8
    log: "{basedir}/pcangsd/global/{file}.log"
    conda:
        "pcangsd"

    shell:
        '''
        conda env list > {log}
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp_weights \
        --sites_save \
        --minMaf 0.05 \
        --threads {threads} \
        --out {params.outdir}/{wildcards.file} \
        &>> {log}
        '''
        
rule run_pcangsd_local:
    input: "{basedir}/angsd/get_maf/{file}.beagle.gz",
    output: 
        done = touch("{basedir}/pcangsd/local/{file}.done"),
    params:
        outdir = "{basedir}/pcangsd/local",
    threads: 8
    log: "{basedir}/pcangsd/local/{file}.log"
    conda:
        "pcangsd"
    shell:
        '''
        mkdir -p {params.outdir}
        ## run pcangsd for PCA only (other options need to be added)
        pcangsd \
        --beagle {input} \
        --snp_weights \
        --sites_save \
        --minMaf 0.05 \
        --threads {threads} \
        --out {params.outdir}/{wildcards.file} \
        &> {log}
        '''