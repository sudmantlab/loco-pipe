import pandas as pd

def extract_bam_column_to_txt(input_tsv_file, output_basedir, separator, filename):
    # Read the TSV file into a DataFrame
    df = pd.read_csv(input_tsv_file, sep=separator)
    
    # Access the 'bam' column
    bam_column = df['bam']
    
    # Write the 'bam' column to a text file
    with open(output_basedir + "/docs/" +  filename, 'w') as txt_file:
        for value in bam_column:
            txt_file.write(str(value) + '\n')

def extract_each_bam_column_to_txt(input_tsv_file, separator, Groupings, base_dir, filename):
    df = pd.read_csv(input_tsv_file, sep=separator)
    for i in pd.unique(pd.Series(df[Groupings])):
        new_df = df[df[Groupings]==i]["bam"]
        with open(base_dir + "/docs/" + i + "_" + filename, "w") as txt_file:
            for each_row in new_df:
                txt_file.write( str(each_row) + '\n')
                
# This rule is to create a global bamlist from the sample table. (By "global", we mean every sample listed in the sample table is being included) 
rule get_bamlist_global:
  input: 
    tsv_file = SAMPLE_TABLE_PATH
  output: 
    txt_file = "{basedir}/docs/"+BAMLIST
  params:
    filename = BAMLIST,
    basedir = BASEDIR
  run:
    extract_bam_column_to_txt(input.tsv_file, params.basedir, "\t", params.filename)

# This rule is to create a bamlist per population from the sample table. (By "local", we mean that we separate samples by population) 
rule get_bamlist_local:
  input: 
    tsv_file = SAMPLE_TABLE_PATH,
  params:
    basedir = BASEDIR,
    grouping = POP_L1_COLNAME,
    filename = BAMLIST
  output: 
    txt_file = expand("{{basedir}}/docs/{population}_"+ BAMLIST, population = POP_L1)
  run:
    extract_each_bam_column_to_txt(input.tsv_file, "\t", params.grouping, params.basedir, params.filename)

# This rule uses samtools inside a conda environment to index the all the input bam files. 
rule index_bam:
  input: 
    bam = "{bam}"
  output: 
    bai = "{bam}.bai"
  threads: 1
  conda:
    "../envs/samtools.yaml"
  shell: "samtools index {input.bam}"
  
# This rule extract the first column of the chromosome table to form a chromosome list. The resulting chromosome list is used to 
# restrict analyses to a predetermined set of chromosomes for certain analyses (e.g. heterozygosity estimation)
rule get_chr_list:
  input: 
      chr_table = "{basedir}/docs/" + config["global"]["chr_table"],
  output:
      chr_list = "{basedir}/docs/chr_list.txt",
  shell:
      '''
      cut -f 1 {input.chr_table} > {output.chr_list}
      '''
