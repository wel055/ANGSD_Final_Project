## ============================================================
## From Raw Read Counts to Relative Expression Strength Measures
## Adapted from 07_practical for GSE95587 Alzheimer's Disease project
## ============================================================

library(ggplot2); theme_set(theme_bw(base_size = 16))
library(magrittr)
library(DESeq2)

## ---- featureCounts results ----

## reading in featureCounts output
df_counts <- read.table("counts/gene_counts.txt", header = TRUE, comment.char = "#")
str(df_counts)

## Clean up sample names: extract SRR IDs from BAM file paths
orig_names <- names(df_counts)
names(df_counts) <- gsub(".*\\.(SRR[0-9]+)\\.Aligned.*", "\\1", orig_names)
## first 6 columns (Geneid, Chr, Start, End, Strand, Length) won't match the regex, kept as-is
str(df_counts)

## ---- Sample mapping ----
## Map SRR IDs to meaningful condition_age labels
sample_info <- data.frame(
  srr = c("SRR5305524", "SRR5305480", "SRR5305508", "SRR5305483", "SRR5305576",
           "SRR5305546", "SRR5305567", "SRR5305591", "SRR5305558", "SRR5305477"),
  condition = c("Control", "Control", "Control", "Control", "Control",
                "AD", "AD", "AD", "AD", "AD"),
  age = c(82, 81, 82, 73, 97,
           82, 82, 81, 73, 97),
  stringsAsFactors = FALSE
)
sample_info$label <- paste0(sample_info$condition, "_", sample_info$age, "_", sample_info$srr)

## ---- countData ----

## gene IDs should be stored as row.names
row.names(df_counts) <- make.names(df_counts$Geneid)

## exclude the columns without read counts (columns 1 to 6)
cts <- df_counts[, -c(1:6)]
head(cts)

## ---- colData ----

## Build colData data.frame matching column order of cts
df_coldata <- data.frame(
  condition = sample_info$condition[match(names(cts), sample_info$srr)],
  age = sample_info$age[match(names(cts), sample_info$srr)],
  row.names = names(cts)
)
df_coldata$condition <- factor(df_coldata$condition, levels = c("Control", "AD"))
df_coldata

## ---- Generate the DESeqDataSet ----

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = df_coldata,
                              design = ~ condition)
dds

## ---- Adding rowData ----
df_rowdata <- df_counts[, 1:6]
rowData(dds) <- df_rowdata
dds

## ---- Accessing counts ----
head(counts(dds))

## Library sizes
colSums(counts(dds))

pdf("counts/library_sizes.pdf", width = 10, height = 6)
barplot(colSums(counts(dds)), las = 2, main = "Library Sizes",
        ylab = "Total read counts", col = ifelse(df_coldata$condition == "AD", "salmon", "lightblue"))
legend("topright", legend = c("AD", "Control"), fill = c("salmon", "lightblue"))
dev.off()

## ---- Remove genes with no reads ----
dim(dds)
keep_genes <- rowSums(counts(dds)) > 0
dds <- dds[keep_genes, ]
dim(dds)

## ---- Normalizing for sequencing depth ----
dds <- estimateSizeFactors(dds)

pdf("counts/sizefactors_vs_libsize.pdf", width = 8, height = 6)
plot(sizeFactors(dds), colSums(counts(dds)),
     ylab = "library sizes", xlab = "size factors", cex = 1.2,
     pch = 16, col = ifelse(colData(dds)$condition == "AD", "salmon", "lightblue"),
     main = "Size Factors vs Library Sizes")
legend("topleft", legend = c("AD", "Control"), pch = 16,
       col = c("salmon", "lightblue"))
text(sizeFactors(dds), colSums(counts(dds)), labels = rownames(colData(dds)),
     pos = 3, cex = 0.7)
dev.off()

## ---- Boxplots: raw vs normalized ----
counts_sf_normalized <- counts(dds, normalized = TRUE)

pdf("counts/boxplots_normalization.pdf", width = 14, height = 6)
par(mfrow = c(1, 2))

boxplot(log2(counts(dds) + 1), notch = TRUE,
        main = "Non-normalized read counts\n(log2)",
        ylab = "log2(read counts)", las = 2, cex = 0.6)

boxplot(log2(counts(dds, normalized = TRUE) + 1), notch = TRUE,
        main = "Size-factor-normalized read counts\n(log2)",
        ylab = "log2(read counts)", las = 2, cex = 0.6)
dev.off()

## ---- Store log-transformed counts ----
assay(dds, "log_counts") <- log2(counts(dds, normalized = FALSE) + 1)
assay(dds, "log_norm_counts") <- log2(counts(dds, normalized = TRUE) + 1)

## ---- Scatterplot of replicates ----
pdf("counts/replicate_scatter.pdf", width = 14, height = 6)
par(mfrow = c(1, 2))

## Pick two Control samples
ctrl_samples <- rownames(df_coldata)[df_coldata$condition == "Control"][1:2]
dds[, ctrl_samples] %>%
  assay(., "log_norm_counts") %>%
  plot(., cex = 0.1, main = paste(ctrl_samples[1], "vs.", ctrl_samples[2], "(Control)"))

## Pick two AD samples
ad_samples <- rownames(df_coldata)[df_coldata$condition == "AD"][1:2]
dds[, ad_samples] %>%
  assay(., "log_norm_counts") %>%
  plot(., cex = 0.1, main = paste(ad_samples[1], "vs.", ad_samples[2], "(AD)"))
dev.off()

## ---- Mean-SD plot (log normalized) ----
pdf("counts/meansd_log_norm.pdf", width = 8, height = 6)
msd_plot <- vsn::meanSdPlot(assay(dds, "log_norm_counts"),
                             ranks = FALSE, plot = FALSE)
print(msd_plot$gg +
  labs(title = "Sequencing depth normalized log2(read counts)",
       y = "standard deviation"))
dev.off()

## ---- rlog transformation ----
dst_rlog <- rlog(dds, blind = TRUE)

## Visual check of rlog
pdf("counts/rlog_comparison.pdf", width = 14, height = 6)
par(mfrow = c(1, 2))
plot(assay(dds, "log_norm_counts")[, 1:2], cex = 0.1,
     main = "size factor and log2-transformed")
plot(assay(dst_rlog)[, 1:2], cex = 0.1,
     main = "rlog transformed",
     xlab = colnames(assay(dst_rlog[, 1])),
     ylab = colnames(assay(dst_rlog[, 2])))
dev.off()

rlog_norm_counts <- assay(dst_rlog)

## ---- Mean-SD plot (rlog) ----
pdf("counts/meansd_rlog.pdf", width = 8, height = 6)
msd_plot <- vsn::meanSdPlot(assay(dst_rlog), ranks = FALSE, plot = FALSE)
print(msd_plot$gg +
  labs(title = "Following rlog transformation",
       x = "Mean", y = "Standard deviation") +
  coord_cartesian(ylim = c(0, 3)))
dev.off()

## ---- Save workspace ----
save.image(file = "counts/AD_normalization.RData")

message("All done! Check counts/ for PDF plots and RData file.")
