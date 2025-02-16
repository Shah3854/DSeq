library(limma)
library(DESeq2)
library(ggplot2)
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(geneplotter)
library(Rsubread)
library(dplyr)
library(genefilter)
library(stringr)
library(LSD)
library(apeglm)
library(tidyverse)
library(matrixStats)
# Load required libraries
library(EnhancedVolcano)

args <- commandArgs(trailingOnly = TRUE)

# Load the count data
counts <- read.csv("Path/PCA_raw.txt", sep = "\t", header = TRUE, row.names = 1)

# Load the condition data from a CSV file
condition_data <- read.csv("Path/condition.csv", header = TRUE)
condition <- as.character(condition_data$Condition)
coldata <- data.frame(row.names = condition_data$Sample, condition)

# Subset the counts data to match the samples in condition data
x <- counts[, rownames(coldata)]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = x[rowSums(x) > 0,], colData = coldata, design = ~condition)

# Run DESeq
ddsDE <- DESeq(dds)

# Get results and filter
results <- results(ddsDE)
res <- na.exclude(as.data.frame(results))
write.csv(results, "path/result.csv")

filter <- res[(abs(res$log2FoldChange) > 1.5 & res$pvalue <= 0.05),]
write.csv(filter, "Path/filter.csv", quote = FALSE, col.names = NA)

# Export normalized counts
normcounts <- counts(ddsDE, normalized = TRUE)
write.csv(normcounts, "Path/norm_counts.csv")

# MA Plot
pdf("Path/maplot.pdf", height = 10, width = 10)
plotMA(results)
dev.off()

# Boxplot
pdf("Path/boxplot.pdf", height = 10, width = 10)
colors <- factor(condition)
boxplot(log2(counts(dds, normalized=FALSE) + 1), col=colors, outline = FALSE,
        main="Box-plot of Normalized Counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample Type", legend=levels(colors), fill=rainbow(length(levels(colors))), cex=0.8)
dev.off()

# Scatter Plot
res <- res %>% mutate(positivecontrol_avg = rowMeans(normcounts[, condition_data$Sample[condition_data$Condition == "PC"]]),
                      treatment_avg = rowMeans(normcounts[, condition_data$Sample[condition_data$Condition == "T"]]))
 res$status <- "NotSignificant"
 res$status[res$log2FoldChange > 0.5] <- "Upregulated"
 res$status[res$log2FoldChange < -0.5] <- "Downregulated"

pdf("Path/scatterplot.pdf", height = 10, width = 10)
 ggplot(res, aes(x = positivecontrol_avg, y = treatment_avg, color = status)) +
   geom_point(alpha = 0.8) +
   scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red", "NotSignificant" = "grey")) +
   labs(title = "Scatter Plot of Normalized Counts (Positive Control vs Treatment)",
        x = "Average Normalized Counts (Positive Control)",
        y = "Average Normalized Counts (Treatment)",
        color = "Gene Status") +
   coord_cartesian(xlim = c(0, 500), ylim = c(0, 500)) +
   theme_minimal()
 dev.off()

# Volcano Plot
pdf("Path/volcano.pdf", height = 10, width = 10)
plot(res$log2FoldChange, -log10(res$pvalue), pch = 20, col = "grey", xlim = c(-10, 10),
     xlab = "log2(FoldChange)", ylab = "-log10(pvalue)", main = "Positive Control vs Treatment")
with(subset(res, pvalue <= 0.05 & log2FoldChange > 2), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))
with(subset(res, pvalue <= 0.05 & log2FoldChange < -2), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
abline(h = -log10(0.05), lty = 2)
abline(v = -2, lty = 2)
abline(v = 2, lty = 2)
dev.off()

# Heatmap
df <- read.csv("Path/norm_counts.csv", row.names=1)
df[df == 0] <- NA
df2 <- df[complete.cases(df),]

top <- df2[order(apply(df2, 1, max), decreasing = TRUE)[1:50],]
write.csv(top, "Path/top50.csv")

heatmap_data <- as.matrix(top)
pdf("Path/heatmap.pdf", height = 10, width = 10)
pheatmap(heatmap_data, cluster_rows=TRUE, cluster_cols=FALSE, show_rownames=TRUE,
         annotation_col=data.frame(Group=condition, row.names=rownames(coldata)),
         color=colorRampPalette(c("red", "white", "purple"))(25))
dev.off()

# Enhanced Volcano Plot
pdf("Path/enhanced_volcano.pdf", height = 10, width = 10)

EnhancedVolcano(res,
                lab = rownames(res),  # Labels for genes
                x = 'log2FoldChange',  # X-axis: log2 fold change
                y = 'pvalue',  # Y-axis: p-value
                title = 'Enhanced Volcano Plot: PC vs Treatment',
                subtitle = 'Differential Expression Analysis',
                xlab = bquote(~Log[2]~ 'Fold Change'),
                ylab = bquote(~-Log[10]~ 'p-value'),
                pCutoff = 0.05,  # Significance threshold
                FCcutoff = 1.5,  # Fold-change cutoff
                pointSize = 2.0,  # Size of points
                labSize = 4.0,  # Label text size
                colAlpha = 0.75,  # Transparency
                legendLabels = c('NS', 'Log2FC', 'p-value', 'Both'),
                legendPosition = 'right',
                col = c('grey30', 'royalblue', 'red2', 'purple'),
                drawConnectors = TRUE,  # Draw lines to labels
                widthConnectors = 0.5
)

dev.off()
