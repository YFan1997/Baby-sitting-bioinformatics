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

### 3. real paper
I choose this paper for recreating protocol: Ochocka, N., Segit, P., Walentynowicz, K.A. et al. Single-cell RNA sequencing reveals functional heterogeneity of glioma-associated brain macrophages. Nat Commun 12, 1151 (2021). https://doi.org/10.1038/s41467-021-21407-w
this is a very comprehensive article, covered the core concepts of single cell analysis : differential expression in different cell type. 
#### 3.1 research goal and result
using scRNA-seq to profile the composition and functions of macriophages(GAMs) in GL261 gilomas in male and female mice, identified the different transcript program in microglia (MG), monocytes/macrophages (Mo/MΦ), and central nervous system (CNS) boarder-associated macrophages (BAMs)
> GL261 is common glioma model
#### 3.2 methods
control: CD11b+ from male and female mice naive brain 


target: CD11b+ from male and female mice tumor brain 


annotation cell cluster


further analyze sub cell type 


identify differential expressed gene, pathway


