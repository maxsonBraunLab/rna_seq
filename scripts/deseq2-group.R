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

if (!dir.exists(snakemake@output$outdir)) {dir.create(snakemake@output$outdir)}

# import data -------------------------------------------------------------------------------------\

# import counts
counts = read.table(snakemake@input$counts, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
gene_names = counts |> pull(Genes)

# import metdata
md = read.table(snakemake@input$md, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
rownames(md) = md$SampleID

# subset md by config/GROUPS.txt if provided
# allow user to upload a contrast combination file for easier interpretation of the output
if (file.exists(snakemake@config$GROUPS)) {

	# import contrast combination file
	groups = read.table(snakemake@config$GROUPS) |> pull(1)

	# subset conditions in md with the conditions in config/GROUPS.txt
	md = md |> filter(Condition %in% groups)
}

# make the md reflect the available sample names
md = md[md$SampleID %in% colnames(counts),]

# make sure md rownames and counts colnames in the same order.
rownames(md) = md$SampleID
counts = counts[, rownames(md)]
stopifnot(rownames(md) == colnames(counts))

# DESeq2 ------------------------------------------------------------------------------------------

# differential expression
dds = DESeqDataSetFromMatrix(countData = counts, colData = md, design = as.formula("~Condition"))
dds.lrt <- DESeq(dds, test="LRT", reduced=~1, parallel = parallel)
res.lrt <- results(dds.lrt, cooksCutoff = Inf, independentFiltering=FALSE)

# Obtain normalized counts
# vsd = vst(dds.lrt, blind = FALSE)
rld <- rlog(dds.lrt, blind=FALSE)

# visualization -----------------------------------------------------------------------------------

colors_fwd <- viridis(100)
colors_rev <- viridis(100, direction = -1)

### sample-sample distances
sampleDists = dist(t(assay(rld)), diag = TRUE)

sampleDistsMatrix = sampleDists |> as.matrix()
rownames(sampleDistsMatrix) <- rld$SampleID
colnames(sampleDistsMatrix) <- NULL

distance_plot = paste0(snakemake@output$outdir, "/", snakemake@config$PROJECT_ID, "-group-sample-dists.pdf")
pdf(distance_plot, width = 12, height = 12)
pheatmap(sampleDistsMatrix,
	fontsize = 12,
	clustering_distance_rows = sampleDists,
	clustering_distance_cols = sampleDists,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	col = colors_rev
)
dev.off()

### PCA plot
pca_plot = paste0(snakemake@output$outdir, "/", snakemake@config$PROJECT_ID, "-group-pca.pdf")
pdf(pca_plot, width = 12, height = 12)
plotPCA(rld, intgroup=c("Condition"))
dev.off()

### Gene expression - top 50 genes post DESeq2 normalization
# top_ge = assay(rld)[top_genes,]
top_gene_indices <- head(order(res.lrt$padj), 50)
top_ge = counts(dds.lrt, normalized=TRUE)[top_gene_indices, ]
rownames(top_ge) = gene_names[top_gene_indices]

top_50 = paste0(snakemake@output$outdir, "/", snakemake@config$PROJECT_ID, "-top-50-genes-heatmap.pdf")
pdf(top_50, width = 12, height = 12)
pheatmap(top_ge,
	main = "Top 50 differential gene expression (DESeq2-normalized counts)",
	scale = "row",
	fontsize = 12,
	cluster_rows = TRUE,
	cluster_cols = TRUE
)
dev.off()

### Gene expression - all genes by individual replicates
all_ge = counts(dds.lrt, normalized=TRUE)
rownames(all_ge) = gene_names
all_ge = all_ge |>
	as.data.frame() |>
	distinct() |>
	as.matrix()
all_ge = all_ge[apply(all_ge, 1, function(x) sd(x) != 0), ]

all_genes_heatmap_individual = paste0(snakemake@output$outdir, "/", snakemake@config$PROJECT_ID, "-all-genes-heatmap-individual.pdf")
pdf(all_genes_heatmap_individual, width = 12, height = 12)
pheatmap(all_ge,
	main = "All gene expression (DESeq2-normalized counts) by individual samples",
	scale = "row",
	fontsize = 12,
	cluster_rows = TRUE,
	cluster_cols = TRUE,
	show_rownames = FALSE
)
dev.off()

### Gene expression - all genes by merged replicates

if (snakemake@config$MERGE_REPLICATES) {

	message("Merging replicates...")

	# read YAML file
	merge_scheme = read_yaml(snakemake@config$REPLICATES)
	for (name in names(merge_scheme)) {

		reps = merge_scheme[[name]]

		if (snakemake@config$MERGE_METHOD == "mean") {
			merged_counts = all_ge[, reps] |> apply(MARGIN = 1, FUN = mean) # same as rowMeans(counts[, reps])
		} else if (snakemake@config$MERGE_METHOD == "median") {
			merged_counts = all_ge[, reps] |> apply(MARGIN = 1, FUN = median)
		} else {
			stop("ERROR: MERGE_METHOD must be either 'mean' or 'median'")
		}

		all_ge = all_ge |>
			as.data.frame() |>
			mutate(!!enexpr(name) := merged_counts) |>
			select(!all_of(reps)) |>
			as.matrix()
		}

	# remove rows with sd=0 one more time.
	all_ge = all_ge[apply(all_ge, 1, function(x) sd(x) != 0), ]

	# define output
	all_genes_heatmap_merged = paste0(snakemake@output$outdir, "/", snakemake@config$PROJECT_ID, "-all-genes-heatmap-merged.pdf")
	pdf(all_genes_heatmap_merged, width = 12, height = 12)
	pheatmap(all_ge,
		main = "All gene expression (DESeq2-normalized counts) by merged samples",
		scale = "row",
		fontsize = 12,
		cluster_rows = TRUE,
		cluster_cols = FALSE, # this respects the order of replicates.yaml
		show_rownames = FALSE
	)
	dev.off()

}