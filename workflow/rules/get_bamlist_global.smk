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
            
# This rule is to create a global bamlist from the metadata table. (By "global", we mean every sample listed in the metadata table is being included) 
rule get_bamlist_global:
  input: 
    tsv_file = META_DATA
  output: 
    txt_file = "{basedir}/docs/"+BAMLIST
  params:
    filename = BAMLIST,
    basedir = BASEDIR
  run:
    extract_bam_column_to_txt(input.tsv_file, params.basedir, "\t",params.filename)
    
  
    