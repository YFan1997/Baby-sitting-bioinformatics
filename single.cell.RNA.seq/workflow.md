## single.cell RNA seq baby sitting
there is an excellent book for single cell RNA-seq:https://www.sc-best-practices.org/preamble.html
we will combine with Seurat R package and take a close look at it
### 1. data structure
scRNA seq is quite similar to the bulk RNA-seq, the key is getting the count from fastq file, so it will also go through alignment and annotation.  Howevere, based on different technique, the count matrix generated from scRNA-seq is different, the unique barcodes would be colnames, while the features would be the rownames, we could use cellranger to process fastq file
```bash
cellranger count --id=sample1 --transcriptome=/path/to/reference --fastqs=/path/to/fastqs --sample=sample1
```

### 2. basic workflow
I recommand visit this website
https://satijalab.org/seurat/articles/pbmc3k_tutorial, which include the basic outline for scRNA data analyze.
briefly, including:

>normalize


>findfeatures


>scaling


>linear dimensional reduction


>clustering the cells


>run non-linear dimensional reduction

