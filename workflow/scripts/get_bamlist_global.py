import pandas as pd

# Read the TSV file into a DataFrame
df = pd.read_csv(snakemake.input.tsv_file, sep='\t')

# Access the 'bam' column
bam_column = df['bam']

# Write the 'bam' column to the output text file
with open(snakemake.output.txt_file, 'w') as txt_file:
    for value in bam_column:
        txt_file.write(str(value) + '\n')
