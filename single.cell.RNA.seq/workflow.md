## single.cell RNA seq baby sitting
there is an excellent book for single cell RNA-seq:https://www.sc-best-practices.org/preamble.html
we will combine with Seurat R package and take a close look at it
### 1. data structure
scRNA seq is quite similar to the bulk RNA-seq, the key is getting the count from fastq file, so it will also go through alignment and annotation.  However, based on different technique, the count matrix generated from scRNA-seq is different, the unique barcodes would be colnames, while the features would be the rownames, we could use cellranger to process fastq file
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

#### 3.3 data preprocessing

scRNA-seq raw data is very large, in this protocol, I will only use two samples and guide through alignment, cellranger. then I will use the provided processed file for further analysis.

```bash
# SRR10004168 :female_control
## downloading may take over 60 minutes, we can prepare the cellranger environment first ( in alignment section)
prefetch SRR10004168 
sam-dump SRR10004168 > SRR10004168.sam
samtools view -bS SRR10004168.sam > SRR10004168.bam
bedtools bamtofastq -i SRR10004168.bam -fq SRR10004168_1.fastq -fq2 SRR10004168_2.fastq
mv SRR10004168_1.fastq female_control_S1_L001_R1_001.fastq
mv SRR10004168_2.fastq female_control_S1_L001_R2_001.fastq
gzip female_control_S1_L001_R1_001.fastq
gzip female_control_S1_L001_R2_001.fastq
```
>important note
in this geo accession, author stored the data in bam format, so fasterq-dump is not suit here, instead, we need to use sam-dump to convert to sam, then bam, then fasterq

>bcl2fastq is normally used for converting single cell to fastq format, we are directly downloading, so it takes different way of care

##### quality check
we don't really need this step, just for remind do quality check every time for fastq file
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

## usecellranger aggr to combine these counts

```
#### download preprocessed data
```bash
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136001/suppl/GSE136001_RAW.tar
mkdir -p processed
tar -xvf GSE136001_RAW.tar -C processed/

# understand file
## here is the cell barcodes which passed quality filtering, this is used to uniquely identify individual cells in the sequencing data
less GSM4039241_f-ctrl-1-filtered-barcodes.tsv.gz | head

AAACCTGAGAGGTTAT-1
AAACCTGCAATGGATA-1
AAACCTGCACAGATTC-1
AAACGGGAGGTGACCA-1
AAACGGGCATCGGACC-1
AAACGGGTCAGGATCT-1
AAACGGGTCAGTTTGG-1
AAACGGGTCGAGGTAG-1
AAACGGGTCTACCTGC-1
AAAGATGCAAGTCTGT-1

## ensmble id, gene symbols
less GSM4039241_f-ctrl-1-filtered-features.tsv.gz | head

ENSMUSG00000051951      Xkr4    Gene Expression
ENSMUSG00000089699      Gm1992  Gene Expression
ENSMUSG00000102343      Gm37381 Gene Expression
ENSMUSG00000025900      Rp1     Gene Expression
ENSMUSG00000025902      Sox17   Gene Expression
ENSMUSG00000104328      Gm37323 Gene Expression
ENSMUSG00000033845      Mrpl15  Gene Expression
ENSMUSG00000025903      Lypla1  Gene Expression
ENSMUSG00000104217      Gm37988 Gene Expression
ENSMUSG00000033813      Tcea1   Gene Expression

## matrix for gene expression counts for each cell
## the 31053 5223 5454144 means there are 31052 rows, each row stands for a gene, 5223 means there are 5223 columns, which there are 5223 cells, while the the 545144 is total number of non-zero elements

## s0976 1 11 means gene at row 309676 has a count of 11 in cell 1

less GSM4039241_f-ctrl-1-filtered-matrix.mtx.gz | head

%%MatrixMarket matrix coordinate integer general
%metadata_json: {"format_version": 2, "software_version": "3.0.1"}
31053 5223 5454144
30976 1 11
30973 1 3
30970 1 7
30969 1 5
30967 1 1
30966 1 7
30964 1 7


## combine data
## this code will make several directory contains barcodes, features, matrix
for sample in GSM4039241_f-ctrl-1 GSM4039242_f-ctrl-2 GSM4039243_f-tumor-1 GSM4039244_f-tumor-2 GSM4039245_m-ctrl-1 GSM4039246_m-ctrl-2 GSM4039247_m-tumor-1 GSM4039248_m-tumor-2; do
    mkdir $sample
    mv ${sample}-filtered-barcodes.tsv.gz $sample/barcodes.tsv.gz
    mv ${sample}-filtered-features.tsv.gz $sample/features.tsv.gz
    mv ${sample}-filtered-matrix.mtx.gz $sample/matrix.mtx.gz
done

```
##### further in R
when get prepared, the major package we will use is Seurat:https://satijalab.org/seurat/ a powerful analyzsis tool for scRNA-seq
