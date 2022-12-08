sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(DESeq2)
library(ggplot2)
library(tibble)
library(dplyr)

parallel <- FALSE
if (snakemake@threads > 1) {
	library("BiocParallel")
	register(MulticoreParam(snakemake@threads))
	parallel <- TRUE
}

# import data --------------------------------------------------------------------------------------

# import counts
counts = read.table(snakemake@input$counts, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# import metdata
md = read.table(snakemake@input$md, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(md) = md$SampleID

# make the md reflect the available sample names
md = md[md$SampleID %in% colnames(counts),]

# make sure md rownames and counts colnames in the same order.
rownames(md) = md$SampleID
counts = counts[, rownames(md)]
stopifnot(rownames(md) == colnames(counts))

# create DESeq2 objects ---------------------------------------------------------------------------

ddsMat = DESeqDataSetFromMatrix(countData = counts, colData = md, design = as.formula("~Condition"))
dds = DESeq(ddsMat, parallel = parallel)

# pairwise analysis -------------------------------------------------------------------------------

# allow user to upload a contrast combination file for easier interpretation of the output
if (!file.exists(snakemake@config$CONTRAST_COMBINATIONS)) {
	contrast_combinations = data.frame(t(combn(unique(md$Condition), 2)))
	colnames(contrast_combinations) = c("c1", "c2")
} else {
	# import contrast combination file
	contrast_combinations = read.table(snakemake@config$CONTRAST_COMBINATIONS, sep = "\t", col.names = c("c1", "c2"))

	# check if all entries are in md
	cc_file_checker = apply( contrast_combinations, 2, function(x) x %in% unique(md$Condition)) |> all()
	if (!cc_file_checker) {
		stop("ERROR: contrast combination file contains a condition not found in the config/metadata.txt file")
	}
}

print("contrast combinations")
print(contrast_combinations)

outdir = snakemake@output$outdir
if (!dir.exists(outdir)) {dir.create(outdir)}

message('writing contrast results...')
for (i in 1:nrow(contrast_combinations)) {
  
  # define contrasts
  c1 = contrast_combinations[i, "c1"]
  c2 = contrast_combinations[i, "c2"]
  contrast_name = paste0(c1, "-vs-", c2)

  # define outputs
  all_file = paste0(outdir, "/", contrast_name, "-all.txt") # contain all gene results for this pw comparison
  sig_file = paste0(outdir, "/", contrast_name, "-", snakemake@params$padj, ".txt") # contain significant gene results
  up_sig = paste0(outdir, "/", contrast_name, "-up-", snakemake@params$padj, ".txt") # contain significantly upregulated genes
  down_sig = paste0(outdir, "/", contrast_name, "-dn-", snakemake@params$padj, ".txt") # contain significantly downregulated genes
  pw_pca = paste0(outdir, "/", contrast_name, "-pca.pdf") # pca of a subset of samples

  # extract contrasts from DESeq object, identify up / down results and intervals.
  temp_res = results(dds, contrast = c("Condition", c1, c2), cooksCutoff = FALSE)
  temp_res_df = temp_res |> 
    as.data.frame() |> 
    rownames_to_column("Genes") |>
    arrange(padj)

  all_sig = temp_res_df |> filter(padj <= snakemake@params$padj)
  upregulated_genes = temp_res_df |> filter(log2FoldChange > 0 & padj <= snakemake@params$padj)
  downregulated_genes = temp_res_df |> filter(log2FoldChange < 0 & padj <= snakemake@params$padj)

  # output PCA plot
  pdf(pw_pca, width = 12, height = 12)
  plotPCA(dds, intgroup=c("Condition"))
  dev.off()

  write.table(all_sig, sig_file, quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(upregulated_genes, up_sig, quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(downregulated_genes, down_sig, quote = FALSE, sep = "\t", row.names = FALSE)

}