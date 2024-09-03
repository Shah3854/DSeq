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
library(dplyr)
library(tidyverse)
library(matrixStats)
library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

# Load the data
counts <- read.csv("/Counts/lnc_counts.tsv" , sep = "\t", header = TRUE, row.names = c(1))
# head(counts)

x<- counts[, c("B_1","B_2", "A_1", "A_2")]
# view(x)

# Define conditions
condition <- c("PC", "PC", "T", "T")
coldata <- data.frame(row.names = colnames(x), condition)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = x[rowSums(x) > 0,], colData = coldata, design = ~condition)

# Run DESeq
ddsDE <- DESeq(dds)

# Get results and filter
results <- results(ddsDE)
res <- na.exclude(as.data.frame(results))

write.csv(results, "/result.csv")

filter <- res[(abs(res$log2FoldChange) > 1.5 & res$pvalue <= 0.05),]
write.csv(filter, "filter.csv", quote = FALSE, col.names = NA)


# Export normalized counts
normcounts <- counts(ddsDE, normalized = TRUE)
write.csv(normcounts, "norm_counts.csv")

#MA PLOT
pdf("maplot.pdf", height = 10, width = 10)
plotMA(results)
dev.off()

#Boxplot
pdf("/boxplot.pdf", height = 10, width = 10)
colors = c(rep("dodgerblue",2), rep("green",2))
boxplot(log2(counts(dds, normalized=FALSE)+1), col=colors, outline = FALSE, main="Box-plot of Normalized counts", xlab="Samples", ylab="log transformed normalized counts")
legend("topright", inset=0, title="Sample type", c("Positive Control", "Treatment"), fill=c("dodgerblue","green"), cex=0.8)
dev.off()
dds <- estimateSizeFactors(dds)


# Calculate averages for scatter plot
res$B_1 <- normcounts[rownames(res), "B_1"]
res$B_2 <- normcounts[rownames(res), "B_2"]
res$A_1 <- normcounts[rownames(res), "A_1"]
res$A_2 <- normcounts[rownames(res), "A_2"]
res$positivecontrol_avg <- (res$B_1 + res$B_2) / 2
res$treatment_avg <- (res$A_1 + res$A_2) / 2

# Define gene status (up, down, NS)
res$status <- "NotSignificant"
res$status[res$log2FoldChange > 0.5] <- "Upregulated"
res$status[res$log2FoldChange < -0.5] <- "Downregulated"

# Create the scatter plot 
pdf("/scatterplot.pdf", height = 10, width = 10)
ggplot(res, aes(x = positivecontrol_avg, y = treatment_avg, color = status)) +
  geom_point(alpha = 0.8) +  # Increase the alpha to make points less transparent
  scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red", "NotSignificant" = "grey")) +  # Use clear, distinct colors
  labs(title = "Scatter Plot of Normalized Counts (Positive Control vs Treatment)",
       x = "Average Normalized Counts (Positive Control)",
       y = "Average Normalized Counts (Treatment)",
       color = "Gene Status") +
  coord_cartesian(xlim = c(0, 500), ylim = c(0, 500)) +  # Adjust the zoom as needed
  theme_minimal()
dev.off()

#multiecdf
# multiecdf(counts(dds, normalized=TRUE)[,], xlab="Meancounts", xlim=c(0,1000))

##volcano-plot 
pdf("/volcano.pdf", height = 10, width = 10)
plot(res$log2FoldChange, -log10(res$pvalue), pch = 20, main = "Positive Control vs Treatment", col = "grey", xlim = c(-10, 10), xlab = "log2(FoldChange)", ylab = "-log10(pvalue)")
with(subset(res, pvalue <= 0.05 & log2FoldChange > 2), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))
with(subset(res, pvalue <= 0.05 & log2FoldChange < -2), points(log2FoldChange, -log10(pvalue), pch = 20, col = "red"))
abline(h = -log10(0.05), lty = 2)
abline(v = -2, lty = 2)
abline(v = 2, lty = 2)
dev.off()


# Heatmap
df <- read.csv("/norm_counts.csv",row.names=NULL)

class(df)
df[df==0] <- NA
df2<-df[complete.cases(df),]

rownames(df2) <- df2[, 1]
df2 <- df2[, -1]

# pick top 50 rows with highest values
top <- df2[order(apply(df2, 1, max), decreasing = TRUE)[1:50],]
write.csv(top, "/top50.csv")
top %>% select(1:4)  -> heatmap_data
heatmap_data %>% pheatmap()
#png("test.png",width=8,height=8,units="in",res=1500)
heatmap_data %>% pheatmap()
heatmap_data %>% log2() -> heatmap_data_log
heatmap_data_log %>% pheatmap()
#png("test2.png",width=8,height=8,units="in",res=1500)
heatmap_data_log %>% pheatmap()
heatmap_data_log - rowMeans((heatmap_data_log)) -> heatmap_data_meanSubtract
heatmap_data_meanSubtract %>% pheatmap()
#png("test3.png",width=8,height=8,units="in",res=1500)
heatmap_data_meanSubtract %>% pheatmap()
#dev.off()

heatmap_data_meanSubtract/rowSds(as.matrix(heatmap_data_log)) -> heatmap_data_zscores
#png("test4.png",width=8,height=8,units="in",res=1500)
heatmap_data_zscores %>% pheatmap(cluster_rows=TRUE , cluster_cols=F , show_rownames = F, border_color = NA)

pdf("/heatmap.pdf", height = 10, width = 10)
annot_cols <- data.frame(
  Group = c( "Positive Control", "Positive Control", "Treatment","Treatment"),
  row.names = colnames(heatmap_data_zscores)
)

color_palette <- colorRampPalette(c("red", "white", "purple"))(25)

pheatmap(
  heatmap_data_zscores,
  show_rownames = F,
  border_color = NA,
  annotation_col = annot_cols,
  cluster_rows = TRUE,
  cluster_cols = F,
  annotation_names_col = T
)
dev.off()



# sessionInfo()
