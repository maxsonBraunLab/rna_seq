import yaml
import pandas as pd
import os

# identify the annotation provided by the STAR index (user-defined)
# filter_col will be used for filtering columns later on.
if snakemake.config["ANNOTATION_STYLE"] == "Genes":
	filter_col = "Gene name lower"
elif snakemake.config["ANNOTATION_STYLE"] == "UCSC":
	filter_col = "UCSC Stable ID"
elif snakemake.config["ANNOTATION_STYLE"] == "ENSEMBL_ID":
	filter_col = "Gene stable ID"
elif snakemake.config["ANNOTATION_STYLE"] == "REFSEQ":
	filter_col = "RefSeq mRNA ID"
elif snakemake.config["ANNOTATION_STYLE"] == "NCBI":
	filter_col = "NCBI gene (formerly Entrezgene) ID"
else:
	print("ERROR: ANNOTATION_STYLE is undefined. Please select from one of the following options: ['Genes', 'UCSC', 'ENSEMBL_ID', 'REFSEQ', 'NCBI']")
	os.exit()

# annotation file -------------------------------------------------------------

# open annotation file
if snakemake.input["annotation"]:
	annotation = pd.read_csv(snakemake.config["ANNOTATION"], sep = "\t")
else:
	print("ERROR: ANNOTATION key in config.yaml is undefined")
	os.exit()

# edit content of the annotation file

# turn gene type and names to lower-case to ensure string matching downstream.
annotation["Gene type lower"] = annotation["Gene type"].transform(lambda x: str.lower(x))
annotation["Gene name lower"] = annotation["Gene name"].transform(lambda x: str(x).lower())

# import biotypes. depending on their values, get them into a list.

biotypes = snakemake.config["BIOTYPES"]
biotype_type = type(biotypes)

if biotype_type == str:
	biotypes = [ snakemake.config["BIOTYPES"] ]
elif biotype_type == list:
	biotypes = snakemake.config["BIOTYPES"]
else:
	print("ERROR: BIOTYPES must be a single entry <string> or multiple entries <list of strings>")

biotypes = [ i.lower() for i in biotypes ]

# filter annotation file contents by biotype
annotation = annotation[annotation["Gene type lower"].isin(biotypes)]
filter_features = list(annotation[filter_col]) # use these features to filter counts table

# open counts file ------------------------------------------------------------

# import counts
counts = pd.read_csv(snakemake.input["counts"], sep = "\t")

print("PREFILTER: Counts table rows = {}".format(counts.shape[0]))
print("PREFILTER: Counts table columns = {}".format(counts.shape[1]))

# filter features
counts["Gene lower"] = counts["Genes"].transform(lambda x: str.lower(x))
counts = counts[counts["Gene lower"].isin(filter_features)]

print("POSTFILTER: Counts table rows = {}".format(counts.shape[0]))
print("POSTFILTER: Counts table columns = {}".format(counts.shape[1] - 1)) # minus one to adjust for lower-case gene names col

counts = counts.drop("Gene lower", axis = 1)

# export filtered counts
counts.to_csv(str(snakemake.output), sep = "\t", index = False)
