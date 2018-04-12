RNASeq Analysis
================

Excel file with read counts produced by running the **RNASeq\_processing.sh** script is imported into RStudio and analyzed for differentially expressed genes (DEGs). The profile of significant DGEs as well as the whole aneuploid population transciptome is then compared to DNA microarrays published in Gasch et al. (2000) Genomic Expression Programs in the Response of Yeast Cells to Environmental Changes. Mol Biol Cell, 11(12): 4241-4257

Step 1: Install/load packages and import sample information
-----------------------------------------------------------

``` r
library(DESeq2)
library(geneplotter)
library(genefilter)
library(gplots)
library(fdrtool)
library(ggplot2)
library(RColorBrewer)
library(openxlsx)
HTSeq.table <- read.xlsx("HTSeq_info.xlsx")  #HTSeq sample information
kable(HTSeq.table)
```

| SampleName | fileName   | genotype  | selection    |  replicate|
|:-----------|:-----------|:----------|:-------------|----------:|
| Teu        | Teu.txt    | euploid   | tetrad       |          1|
| Taneu1     | Taneu1.txt | aneuploid | tetrad       |          2|
| Taneu2     | Taneu2.txt | aneuploid | tetrad       |          2|
| Aeu        | Aeu.txt    | euploid   | a\_selection |          3|
| Aaneu1     | Aaneu1.txt | aneuploid | a\_selection |          4|
| Aaneu2     | Aaneu2.txt | aneuploid | a\_selection |          4|
| Aaneu3     | Aaneu3.txt | aneuploid | a\_selection |          4|

Step 2: Build DESeq dataset from HTSeq data
-------------------------------------------

``` r
DESeq.data <- DESeqDataSetFromHTSeqCount(sampleTable = HTSeq.table, design = ~genotype + 
    selection)
colData(DESeq.data)$genotype <- factor(colData(DESeq.data)$genotype, levels = c("euploid", 
    "aneuploid"))
colData(DESeq.data)$selection <- factor(colData(DESeq.data)$selection)
```

Step 3: Normalize samples and filter out low-count transcripts
--------------------------------------------------------------

``` r
DESeq2.data <- estimateSizeFactors(DESeq.data)
# Keep only transcripts with counts greater than 5 in altleast 3 samples
counts.filter <- rowSums(counts(DESeq2.data, normalized = TRUE) >= 5) >= 3
DESeq2.data <- DESeq2.data[counts.filter, ]
# Check count distribution after normalization
multiecdf(counts(DESeq2.data, normalized = TRUE), xlim = c(0, 5000), xlab = "Mean Counts")
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/normalize-1.png)

``` r
multidensity(counts(DESeq2.data, normalized = TRUE), xlim = c(0, 5000), xlab = "Mean Counts")
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/normalize-2.png)

Step 4: Assess data quality and sample similarity
-------------------------------------------------

``` r
rld <- rlogTransformation(DESeq2.data, blind = TRUE)  #log-transform data for QC only
# Sample distance heatmap
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, trace = "none", col = rev(hmcol), margin = c(10, 10), Rowv = FALSE, 
    symm = TRUE, density.info = "none", key.title = NA, key.xlab = "Sample Distance", 
    dendrogram = "none")
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/QC-1.png)

``` r
# PCA
DESeq2::plotPCA(rld, intgroup = c("genotype", "selection")) + theme_bw() + scale_color_manual(values = c("cornflowerblue", 
    "indianred", "blue", "red")) + theme(panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), axis.title = element_text(size = 16), 
    legend.position = "right") + guides(color = guide_legend("Strain"))
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/QC-2.png)

Step 5: Differential expression analysis
----------------------------------------

``` r
DESeq2.data <- estimateDispersions(DESeq2.data)
plotDispEsts(DESeq2.data)  #Mean-dispersion relationship
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/DESeq-1.png)

``` r
DESeq2Table <- nbinomWaldTest(DESeq2.data)
DESeq2Res <- results(DESeq2Table, pAdjustMethod = "BH", alpha = 0.05, contrast = c("genotype", 
    "aneuploid", "euploid"))
hist(DESeq2Res$padj, main = "Adjusted p-value distribution", xlab = "DESeq2 p-values", 
    breaks = seq(0, 1, 0.05))
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/DESeq-2.png) \#\# Step 6: Estimate FDR using empirical null distribution Check if the enrichment of low p-values is because the null variance is underestimated (this could give rise to false positives).

``` r
DESeq2Res <- DESeq2Res[!is.na(DESeq2Res$padj), ]  #Remove genes with NA p-values
DESeq2Res <- DESeq2Res[!is.na(DESeq2Res$pvalue), ]
FDR.DESeq2Res <- fdrtool(DESeq2Res$stat, plot = FALSE, verbose = FALSE)
kable(FDR.DESeq2Res$param[, 3:6])  #SD is 2.7524, greater than 1
```

|         |           |
|:--------|----------:|
| eta0    |  1.0000000|
| eta0.SE |  0.0012862|
| sd      |  2.7531433|
| sd.SE   |  0.0325164|

``` r
DESeq2Res <- DESeq2Res[, -which(names(DESeq2Res) == "padj")]  #Remove adjusted p-value column
DESeq2Res[, "padj"] <- FDR.DESeq2Res$pval
hist(DESeq2Res$padj, main = "Adjusted p-value distribution", xlab = "DESeq2 p-values", 
    breaks = seq(0, 1, 0.05))
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/Empirical%20null-1.png)

``` r
RNASeq_FoldChange <- data.frame(DESeq2Res)  #Tabulated log2FoldChange(aneuploid/euploid) 
# Differentially expressed genes (DEGs)
DEGs <- RNASeq_FoldChange[RNASeq_FoldChange$padj < 0.05, ]
DEGs$Gene <- row.names(DEGs)
row.names(DEGs) <- NULL
kable(head(DEGs, 4))
```

|    baseMean|  log2FoldChange|      lfcSE|      stat|  pvalue|       padj| Gene    |
|-----------:|---------------:|----------:|---------:|-------:|----------:|:--------|
|  28606.4667|       1.9485554|  0.2788650|  6.987452|       0|  0.0111491| YAL003W |
|    892.5034|       1.9643507|  0.3424205|  5.736663|       0|  0.0371895| YAL019W |
|  10731.5152|       0.5126114|  0.0882116|  5.811158|       0|  0.0347950| YAL023C |
|    664.2886|       3.8562211|  0.6376223|  6.047814|       0|  0.0280423| YAL025C |

``` r
# MA plot
ggplot() + geom_point(data = RNASeq_FoldChange, aes(x = log10(baseMean), y = log2FoldChange), 
    color = "darkgrey", size = 0.8) + geom_point(data = RNASeq_FoldChange[RNASeq_FoldChange$padj <= 
    0.05, ], aes(x = log10(baseMean), y = log2FoldChange), color = "red2", size = 0.8) + 
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") + scale_x_continuous(breaks = c(1, 
    3, 5), labels = c(expression(10^{
    1
}), expression(10^{
    3
}), expression(10^{
    5
}))) + ylab(expression(~log[2] ~ (Expression ~ Fold ~ Change))) + theme_bw() + 
    xlab("Basal Mean Expression") + theme(panel.grid = element_blank(), axis.text = element_text(size = 16), 
    axis.title = element_text(size = 18))
```

![](RNASeq_DGEAnalysis_files/figure-markdown_github-ascii_identifiers/Empirical%20null-2.png)

Step 7: Correlation with microarrays in Gasch et al. (2000)
-----------------------------------------------------------

``` r
Gasch.data <- read.xlsx("Gasch_data.xlsx")  #Load Gasch data
Gasch.data <- Gasch.data[, -grep("000", colnames(Gasch.data))]  #Remove WT unstressed reference conditions
Gasch.data.DEGs <- merge(DEGs, Gasch.data, by.x = "Gene", by.y = "UID")
Correlation.matrix <- data.frame(Conditions = colnames(Gasch.data.DEGs)[c(3, 
    10:178)], cor(Gasch.data.DEGs[, c(3, 10:178)], use = "pairwise.complete.obs", 
    method = "spearman"))
Correlation.aneuploidy <- Correlation.matrix[order(Correlation.matrix$log2FoldChange, 
    decreasing = TRUE), 1:2]
row.names(Correlation.aneuploidy) <- NULL
kable(head(Correlation.aneuploidy, 7), align = "c")
```

|          Conditions         | log2FoldChange |
|:---------------------------:|:--------------:|
|        log2FoldChange       |    1.0000000   |
| Hypo-osmotic.shock.-.15.min |    0.7199484   |
| Hypo-osmotic.shock.-.30.min |    0.7132712   |
|  37C.to.25C.shock.-.15.min  |    0.6992636   |
|  37C.to.25C.shock.-.30.min  |    0.6828444   |
|  37C.to.25C.shock.-.45.min  |    0.6180003   |
|      dtt.015.min.dtt-2      |    0.6043048   |
