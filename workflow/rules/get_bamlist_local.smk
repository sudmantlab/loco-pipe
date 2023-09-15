import pandas as pd
def extract_each_bam_column_to_txt(input_tsv_file, separator, Groupings, base_dir, filename = "bamlist_min1x.txt"):
    df = pd.read_csv(input_tsv_file, sep=separator)
    for i in pd.unique(pd.Series(df[Groupings])):
        new_df = df[df[Groupings]==i]["bam"]
        with open(base_dir + "/docs/" + i + "_" + filename, "w") as txt_file:
            for each_row in new_df:
                txt_file.write( str(each_row) + '\n')


rule get_bamlist_local:
  input: 
    tsv_file = META_DATA,
  params:
    basedir = BASEDIR,
    grouping = POP_L1,
    filename = BAMLIST
  output: 
    txt_file = "{basedir}/docs/{population}_"+ BAMLIST
  run:
    extract_each_bam_column_to_txt(input.tsv_file, "\t", params.grouping, params.basedir, params.filename)
    