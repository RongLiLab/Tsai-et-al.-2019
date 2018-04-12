DNASeq Copy Number Analysis
================

Excel file with read counts produced by running the **DNASeq\_processing.sh** script is imported into RStudio and analyzed for copy number variants.

Step 1: Install/load packages and import data
---------------------------------------------

``` r
library(DNAcopy)
library(GenomicRanges)
library(cn.mops)
library(openxlsx)
library(ggplot2)

read.counts <- read.xlsx("DNASeq_Counts.xlsx")  #Import data
kable(head(read.counts, n = 4))
```

| chr |  start|   end|  width| strand |  haploid|    Teu|  Taneu1|  Taneu2|    Aeu|  Aaneu1|  Aaneu2|  Aaneu3|
|:----|------:|-----:|------:|:-------|--------:|------:|-------:|-------:|------:|-------:|-------:|-------:|
| I   |      1|  1501|   1501| \*     |     3439|   4711|    3894|    3907|   5506|    6425|    5549|    4044|
| I   |   1502|  2502|   1001| \*     |     2407|   3796|    1657|    2563|   1356|    3463|    2109|    1806|
| I   |   2503|  4003|   1501| \*     |    17696|  27403|   19891|   22985|  26875|   35133|   23248|   14258|
| I   |   4004|  5504|   1501| \*     |    23969|  38736|   25937|   27375|  31836|   41929|   28679|   23596|

Step 2: Normalize read counts
-----------------------------

``` r
reads.gr <- makeGRangesFromDataFrame(read.counts, keep.extra.columns = TRUE)  #Create GRanges object
normalized.counts <- data.frame(normalizeGenome(reads.gr, normType = "median", 
    ploidy = rep(1, 8)))  #Median normalization
kable(head(normalized.counts, n = 4))
```

| seqnames |  start|   end|  width| strand |    haploid|        Teu|     Taneu1|     Taneu2|        Aeu|     Aaneu1|     Aaneu2|     Aaneu3|
|:---------|------:|-----:|------:|:-------|----------:|----------:|----------:|----------:|----------:|----------:|----------:|----------:|
| I        |      1|  1501|   1501| \*     |   3518.063|   4233.738|   4047.019|   3995.528|   5243.486|   5157.068|   5428.717|   4633.443|
| I        |   1502|  2502|   1001| \*     |   2462.337|   3411.435|   1722.114|   2621.075|   1291.349|   2779.600|   2063.284|   2069.238|
| I        |   2503|  4003|   1501| \*     |  18102.832|  24626.859|  20672.639|  23505.816|  25593.659|  28199.733|  22744.062|  16336.208|
| I        |   4004|  5504|   1501| \*     |  24520.048|  34811.737|  26956.224|  27995.288|  30318.130|  33654.587|  28057.337|  27035.290|

Step 3: Transform data to log2 ratio relative to haploid control
----------------------------------------------------------------

``` r
log2ratios <- normalized.counts
log2ratios[, grep("eu", colnames(log2ratios))] <- log2(normalized.counts[, grep("eu", 
    colnames(normalized.counts))]/normalized.counts$haploid)
log2ratios[, 7:13] <- apply(log2ratios[, 7:13], 2, function(x) {
    x[(!is.finite(x))] <- 0  #Convert log2(x/0) and log2(0/0) to 0.
    x
})
log2ratios <- log2ratios[, -6]  #Remove haploid data
kable(head(log2ratios, n = 4))
```

| seqnames |  start|   end|  width| strand |        Teu|      Taneu1|     Taneu2|         Aeu|     Aaneu1|      Aaneu2|      Aaneu3|
|:---------|------:|-----:|------:|:-------|----------:|-----------:|----------:|-----------:|----------:|-----------:|-----------:|
| I        |      1|  1501|   1501| \*     |  0.2671509|   0.2020785|  0.1836051|   0.5757451|  0.5517700|   0.6258299|   0.3973033|
| I        |   1502|  2502|   1001| \*     |  0.4703506|  -0.5158478|  0.0901303|  -0.9311494|  0.1748489|  -0.2550858|  -0.2509289|
| I        |   2503|  4003|   1501| \*     |  0.4440173|   0.1915072|  0.3768023|   0.4995710|  0.6394661|   0.3292746|  -0.1481423|
| I        |   4004|  5504|   1501| \*     |  0.5056120|   0.1366566|  0.1912222|   0.3062189|  0.4568413|   0.1944163|   0.1408820|

Step 4 (a): Copy number analysis using DNAcopy
----------------------------------------------

``` r
CNA.data <- CNA(genomdat = log2ratios[, 6:12], chrom = log2ratios$seqnames, 
    maploc = log2ratios$start, data.type = "logratio", sampleid = colnames(log2ratios)[6:12])

CNA.smoothed <- smooth.CNA(CNA.data)
segments.CNA <- DNAcopy::segment(CNA.data, verbose = 0, min.width = 2)
segments.data <- segments.CNA$output
kable(head(segments.data, n = 4))
```

| ID  | chrom |  loc.start|  loc.end|  num.mark|  seg.mean|
|:----|:------|----------:|--------:|---------:|---------:|
| Teu | I     |          1|   228680|       180|    0.0938|
| Teu | II    |          1|   121099|        99|    0.0016|
| Teu | II    |     121600|   466864|       265|   -0.1514|
| Teu | II    |     468365|   472868|         4|   -0.6969|

Step 4 (b) : Copy number analysis using Mixture of Poissons (cn.mops)
---------------------------------------------------------------------

``` r
copy.number <- haplocn.mops(reads.gr, norm = 1, minWidth = 1, I = c(0.025, 1:8), 
    classes = paste("CN", 0:8, sep = ""), priorImpact = 0.8, parallel = 4)
resCNMOPS <- calcIntegerCopyNumbers(copy.number)  #Integer Copy nmbers of all segments
variants.mops <- as.data.frame(resCNMOPS@integerCopyNumber)
metadata <- strsplit(row.names(variants.mops), "_")
CNA.mops <- data.frame(Chromosome = sapply(metadata, "[", 1), Start = sapply(metadata, 
    "[", 2), End = sapply(metadata, "[", 3), variants.mops, row.names = NULL)
kable(head(CNA.mops, n = 4))
```

| Chromosome | Start | End  | haploid | Teu | Taneu1 | Taneu2 | Aeu | Aaneu1 | Aaneu2 | Aaneu3 |
|:-----------|:------|:-----|:--------|:----|:-------|:-------|:----|:-------|:-------|:-------|
| I          | 1     | 1501 | CN1     | CN1 | CN1    | CN1    | CN1 | CN1    | CN1    | CN1    |
| I          | 1502  | 2502 | CN1     | CN1 | CN1    | CN1    | CN1 | CN1    | CN1    | CN1    |
| I          | 2503  | 4003 | CN1     | CN1 | CN1    | CN1    | CN1 | CN1    | CN1    | CN1    |
| I          | 4004  | 5504 | CN1     | CN1 | CN1    | CN1    | CN1 | CN1    | CN1    | CN1    |
