import pandas as pd


"""Function accepts a STAR output directory and compiles all sample information from ReadsPerGene.out.tab
 
Args:
    snakemake.input (list): list of globbed wildcards STAR ReadsPerGene.out.tab
    project_title (str): Project title for compiled STAR counts
Returns:
    Compiled STAR counts as tab delimited file.
"""

colnames = snakemake.params.samples
tables = [pd.read_csv(fh, header=None, skiprows=4, usecols=[0,3], index_col=0, sep = '\t',  names = ['Genes',fh.split('/')[-2].split('_')[0]]) for fh in snakemake.input]
joined_table = pd.concat(tables, axis=1)
joined_table.columns = colnames
joined_table_sorted = joined_table.reindex(sorted(joined_table.columns), axis = 1)
joined_table_sorted.to_csv(snakemake.output[0], sep='\t')