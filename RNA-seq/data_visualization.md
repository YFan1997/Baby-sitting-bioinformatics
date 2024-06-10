## let's start doing data visulaization by R
### 1. load necessary r packages
```{r}
# BiocManager::install("DESeq2")
# BiocManager::install("pheatmap")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("edgeR")
# BiocManager::install(c("tximport","readr","biomaRt"))
library(DESeq2)
library(edgeR)
library(tximport)
library(pheatmap)
library(biomaRt)
library(EnhancedVolcano)
```
we actually mainly use DESeq2, but other packages are also very useful for deeper understanding and management of our result.
In DESeq2, there are many embeded algorithms which users can directly use, but, I really want to digging in the meaning of these parameters, we will find these would help a lot in understanding plots and make our own plots based on specific needs.
> log2FC
this parameter stands for log2 fold change. firstly, fold change measures how much quantity changes between two conditions.let's say condition A in gene x's expression level(read count) is a, and condition B in gene x' expression level(read count) is b, then the fold change will be b/a, the log2FC will be log2(b/a). Why using log2, as it could make data symmetric, as if b = a, log2(b/a) = 0, etc.  this concept will help us understanding plotMA




