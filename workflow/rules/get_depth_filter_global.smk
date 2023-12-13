# This rule reads in the depthGlobal files, fits them into a normal distribution, 
# establishes min and max filters based on the mean and standard deviation of the fitted distribution.
# These filters will be used by the next steps of the analysis (e.g. get_site_list_global.smk and snp_calling_global.smk)
# This rule outputs a density plot along with a tsv file on read depth distribution and min&max filter.

rule get_depth_filter_global:
    input: expand("{{basedir}}/angsd/get_depth_global/{chr}.depthGlobal", chr = CHRS)
    output: 
        "{basedir}/angsd/get_depth_global/depth_filter.tsv",
        "{basedir}/figures/depth/depth_filter.png"
    threads: 4
    params:
        indir = "{basedir}/angsd/get_depth_global/",
        plot_dir = "{basedir}/figures/depth/",
        chr_list = CHR_TABLE_PATH,
        n_sd = config["get_depth_filter_global"]["n_sd"],
        rscript = config["global"]["scriptdir"] + "/get_depth_filter.R"
    log: "{basedir}/angsd/get_depth_global/get_depth_filter.log"
    conda: "../envs/r.yaml"
    shell:
        '''
        mkdir -p {params.plot_dir}
        Rscript {params.rscript} {params.indir} {params.plot_dir} {params.chr_list} {params.n_sd} &> {log}
        '''




