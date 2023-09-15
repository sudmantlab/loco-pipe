# This rule uses samtools inside a conda environment 
# to index the all the input bam files. 
rule index_bam:
  input: 
    bam = "{bam}"
  output: 
    bai = "{bam}.bai"
  threads: 1
  conda:
    "../envs/samtools.yaml"
  shell: "samtools index {input.bam}"
