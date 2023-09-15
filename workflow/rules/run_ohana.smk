# This snakemake file runs Ohana to conduct admixture analysis to characterize population structure. 
# More detailed description of Ohana could be found here: https://github.com/jade-cheng/ohana/blob/master/README.md.

# This rule converts the subsetted beagle formatted genotype likelihood file to an .lgm format file 
# which is needed for Ohana.
rule get_ohana_input_global:
    input: rules.combine_subsetted_beagle_global.output.subsetted_beagle
    output: "{basedir}/ohana/global/combined.subsetted.lgm"
    params: outdir = "{basedir}/ohana/global"
    conda: "../envs/ohana.yaml"
    threads: 1
    shell:
        '''
        mkdir -p {params.outdir}
        convert bgl2lgm {input} {output}
        '''
        
# After the file conversion from the get_ohana_input_global rule, this rule applies "qpas"
# method in Ohana with different hyperparameters for k. You could adjust the range for k by 
# changing the "min_k" and "max_k" in the config file. The final outputs are q_matrix and f_matrix. 
# The former records the admixture proportions for each sample, and the latter records the allele 
# frequencies in each ancestral population.
# A well-explained example could be found here: https://github.com/jade-cheng/ohana/blob/master/README.md.
rule run_ohana_global:
     input: "{basedir}/ohana/global/combined.subsetted.lgm"
     output: 
         q_mat = "{basedir}/ohana/global/combined.subsetted.k{k}.q.matrix",
         f_mat = "{basedir}/ohana/global/combined.subsetted.k{k}.f.matrix",
         done = touch("{basedir}/ohana/global/combined.subsetted.k{k}.done"),
     threads: 4
     log: "{basedir}/ohana/global/combined.subsetted.k{k}.log"
     conda: "../envs/ohana.yaml"
     shell:
         '''
         qpas {input} \
         -k {wildcards.k} \
         -qo {output.q_mat} \
         -fo {output.f_mat} \
         -mi 50 &> {log}
         '''
# The next two rules inherit the same logic from get_ohana_input_global and run_ohana_global rules,
# but these rules uses data specific to one population.
rule get_ohana_input_local:
    input: rules.combine_subsetted_beagle_local.output.subsetted_beagle
    output: "{basedir}/ohana/local/{population}.combined.subsetted.lgm"
    params: outdir = "{basedir}/ohana/local"
    conda: "../envs/ohana.yaml"
    threads: 1
    shell:
        '''
        mkdir -p {params.outdir}
        convert bgl2lgm {input} {output}
        '''

rule run_ohana_local:
     input: "{basedir}/ohana/local/{population}.combined.subsetted.lgm"
     output: 
         q_mat = "{basedir}/ohana/local/{population}.combined.subsetted.k{k}.q.matrix",
         f_mat = "{basedir}/ohana/local/{population}.combined.subsetted.k{k}.f.matrix",
         done = touch("{basedir}/ohana/local/{population}.combined.subsetted.k{k}.done"),
     threads: 4
     log: "{basedir}/ohana/local/{population}.combined.subsetted.k{k}.log"
     conda: "../envs/ohana.yaml"
     shell:
         '''
         qpas {input} \
         -k {wildcards.k} \
         -qo {output.q_mat} \
         -fo {output.f_mat} \
         -mi 50 &> {log}
         '''
