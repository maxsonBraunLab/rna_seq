sink(file(snakemake@log[[1]], open = "wt"), type = "message")

library(DESeq2)
library(ggplot2)
library(tibble)
library(dplyr)
library(viridisLite)
library(pheatmap)
library(yaml)

parallel <- FALSE
if (snakemake@threads > 1 & snakemake@params$use_singularity == "false") {
	library("BiocParallel")
	register(MulticoreParam(snakemake@threads))
	parallel <- TRUE
}

c1 = snakemake@params$c1
c2 = snakemake@params$c2
outdir = snakemake@params$outdir
contrast_name = paste0(snakemake@params$c1, '-vs-', snakemake@params$c2)

# import data -------------------------------------------------------------------------------------

# import counts
counts = read.table(snakemake@input$counts, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = "Genes")

# import metdata
md = read.table(snakemake@input$md, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# subset and counts matrix by contrast groups
md = md |> filter(Condition %in% c(c1, c2))
rownames(md) = md$SampleID

# make sure md rownames and counts colnames in the same order.
counts = counts[, rownames(md)]
stopifnot(rownames(md) == colnames(counts))

# DESeq2 ------------------------------------------------------------------------------------------

# differential expression
ddsMatrix = DESeqDataSetFromMatrix(countData = counts, colData = md, design = as.formula("~Condition"))
dds <- DESeq(ddsMatrix, parallel = parallel)

# Obtain normalized counts
# vsd = vst(dds.lrt, blind = FALSE)
rld <- rlog(dds, blind=FALSE)

# write outputs -----------------------------------------------------------------------------------

# define outputs
all_file = paste0(outdir, "/", contrast_name, "-all.txt") # contain all gene results for this pw comparison
sig_file = paste0(outdir, "/", contrast_name, "-", snakemake@params$padj, ".txt") # contain significant gene results
pw_pca = paste0(outdir, "/", contrast_name, "-pca.pdf") # pca of a subset of samples
pw_ma = paste0(outdir, "/", contrast_name, "-MA.pdf") # MA of a subset of samples

# extract results from DESeq object
results = results(dds, contrast = c("Condition", c1, c2), cooksCutoff = FALSE)
results_df = results |>
	as.data.frame() |>
	rownames_to_column("Genes") |>
	arrange(padj)
all_sig = results_df |> filter(padj <= snakemake@params$padj)

# export tables
write.table(results_df, all_file, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(all_sig, sig_file, quote = FALSE, sep = "\t", row.names = FALSE)

# output PCA plot
pdf(pw_pca, width = 12, height = 12)
plotPCA(rld, intgroup=c("Condition"))
dev.off()

# plot MA
pdf(pw_ma, width = 12, height = 12)
plotMA(results, intgroup=c("Condition"))
dev.off()
