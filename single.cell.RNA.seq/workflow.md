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
#### 3.2 workflow
>control: CD11b+ from male and female mice naive brain 
>target: CD11b+ from male and female mice tumor brain 


annotation cell cluster


further analyze sub cell type 


identify differential expressed gene, pathway

#### 3.3 data proprocessing

scRNA-seq raw data is very large, in this protocol, I will only use two samples and guide through alignment, cellranger. then I will use the provided processed file for further analysis.

```bash
# SRR10004168 :female_control
# SRR10004170 :female_tumor
## downloading may take over 60 minutes, we can prepare the cellranger environment first ( in alignment section)
prefetch SRR10004168 && fasterq-dump 
prefetch SRR10004170 && fasterq-dump 

mv SRR10004168.fastq female_control.fastq
mv SRR10004170.fastq female_tumor.fastq
## the downloaded is merged version, we can use seqtk to seperate 
seqtk seq -1 female_control.fastq > female_control_R1.fastq
seqtk seq -2 female_control.fastq > female_control_R2.fastq
# do the same for female_tumore

## to match the name convention, using gzip and further rename
gzip female_control_R1.fastq
gzip female_control_R2.fastq
mv female_control_R1.fastq.gz female_control_S1_L001_R1_001.fastq.gz
mv female_control_R2.fastq.gz female_control_S1_L001_R2_001.fastq.gz
## do the similiar for female_tumor

```
>bcl2fastq is normally used for converting single cell to fastq format, we are directly downloading, so it takes different way of care

##### quality check
```bash
fastqc female_control_S1_L001_R1_001.fastq.gz female_control_S1_L001_R2_001.fastq.gz -o fastqc_output
```
##### scRNA-seq alignment
here the major software used is cellranger:https://www.10xgenomics.com/support/software/cell-ranger/latest
downloading could follow this protocol: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in#tutorial
add to working envrionment
and I use mainly cellranger count: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct
ref: mm10 prebuid in cellranger
```bash
export PATH=/mnt/pv_compute/yifan/tools/cellranger/cellranger-9.0.0:$PATH
mm10_ref="/mnt/pv_compute/yifan/tools/cellranger/mm10/refdata-gex-GRCm39-2024-A"
fastq_path="/mnt/pv_compute/yifan/practice/scRNA.practice/raw.data/fastqs"

cellranger count --id=female_control \
                 --transcriptome=$mm10_ref \
                 --fastqs=$fastq_path \
                 --sample=female_control \
                 --localcores=8 \
                 --localmem=32 \
                 --create-bam=true

cellranger count --id=female_tumor \
                 --transcriptome=$mm10_ref \
                 --fastqs=$fastq_path \
                 --sample=female_tumor \
                 --localcores=8 \
                 --localmem=32 \
                 --create-bam=true

## usecellranger aggr to combine these counts

```
