library(DESeq2)

# import data -------------------------------------------------------------------------------------

# import counts
counts = read.table(snakemake@input$counts, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = "Genes")

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
ddsMat = DESeqDataSetFromMatrix(countData = counts, colData = md, design = as.formula(snakemake@config$MODEL))
dds = DESeq(ddsMat)

message('calculating deseq...')

# Export DESeq2-normalized counts -----------------------------------------------------------------
normCounts = counts(dds, normalized=TRUE)
normCounts = 0.1 + normCounts # ensure nonzero
log2normCounts = log2(normCounts)
log2normCounts = as.data.frame(log2normCounts)
write.table(normCounts, snakemake@output$norm_counts, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
write.table(log2normCounts, snakemake@output$log2_norm_counts, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
