# RNA-Seq 
RNA-seq is widely used next-generation sequencing method to determine the association between RNA and biosamples,
I will start from the wet lab process and then introduce a practical RNA-seq analyze procedure.

ref1: Kukurba KR, Montgomery SB. RNA Sequencing and Analysis. 
Cold Spring Harb Protoc. 2015;2015(11):951-969. Published 2015 Apr 13. doi:10.1101/pdb.top084970
ref2: Koch CM, Chiu SF, Akbarpour M, et al. A Beginner's Guide to Analysis of RNA Sequencing Data. 
Am J Respir Cell Mol Biol. 2018;59(2):145-157. doi:10.1165/rcmb.2017-0430TR

## wet lab process
> prepare samples, and extract RNA
>
> 
in this step, researchers often use the ratios of 28s/18s ribosomal bands to estimates sample intergrity.
> creation of an RNA-seq library
>
> 
isolating the target RNA molecule, and convert RNA to cDNA, ligating with sequencing
adapter, there are multiple methods for library design based on different needs.
> quantification


"spike-in" sample was added to normalized the coverage and sensitivity.
### understanding the how to conduct RNA-seq is essential even though it's not our focus in this topic, it would help us trouble shooting from the analyze result.

## RNA seq data analyze
the standard original input analyze file for RNA-seq is fastq, normally, align the sequence to the reference genome and get reads counts is the first step for analyzing fastq data.
Nowdays, there are many open-sourced tool could be used for alignments, when we make a decision, we also need to consider the features of the software.
The following procedure are based on Koch CM, Chiu SF, Akbarpour M, et al

### 1. obtain the data
check the file for instruction of data downloading, and the data we are using is:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116583
click on SRA Run Selector in the bottom of page, and select all, click AssessionList, it would download a txt file which contains the sample name we want.
```bash
# prepare for download
mkdir RNA-seq
mkdir raw_data
cd RNA-seq/raw_data
# here we need to download 12 files, if the server have multiple threads, we can directly use the following command, if we are going to huge more files for processing, add -P8 after parallel to limit the usage of the threads. nohup could make your command hang in the server without termination
# this downloading process would tabke more than 20 minutes
nohup parallel 'prefetch {} && fasterq-dump {}' :::: SRA_download_list &

# while we are waiting for downloading, let's write some command line to rename the fastq file.
# create a file called mapping_list, which contain: use
nano mapping_list
###
R05_AM_Allo24h
R06_AM_Allo24h
R03_AM_Allo2h
R04_AM_Allo2h
R01_AM_Allo2h
R02_AM_Allo2h
N03_AM_Naive
N04_AM_Naive
N01_AM_Naive
N02_AM_Naive
R07_AM_Allo24h
R08_AM_Allo24h
###
```
```bash
# write a script, called rename.sh

#!/bin/bash
#{rename.sh}
# Read the lists into arrays
mapfile -t sra_list < SRA_download_list
mapfile -t name_list < mapping_list

# Loop through each pair of entries and rename the files
for i in "${!sra_list[@]}"; do
    sra_id="${sra_list[$i]}"
    new_name="${name_list[$i]}"
    
    # Rename the fastq file
    mv "${sra_id}.fastq" "${new_name}.fastq"
done

# run this script
bash rename.sh
```
### 2.investigae data
here we mainly use fastqc method
```bash
mkdir fastqc_output
nohup parallel fastqc {}.fastq -o fastqc_output/ :::: mapping_list &
```
got html file, check adapter content, if they are green, ready for next step, if not consider trimming method such as cutadapt,trimmomatic etc.

### 3.alignment
in this experiment, the research are using mm10 for reference gene
```bash
#download and unzip reference genome
wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
# download and unzip annotation
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz
gunzip Mus_musculus.GRCm38.99.gtf.gz

# build index by hisat2
# again, there are many software for doing the similar job 

nohup hisat2-build mm10.fa mm10_index &

#some small step make life easier
# cp mapping_list alignment
# mv mapping_list name
# alignment
cd /mnt/pv_compute/yifan/practice/RNA-seq/raw_data/alignment
nohup parallel "hisat2 -p 8 -x path/to/mm10_index -U path/to/{}.fastq -S path/to/sam/{}.sam" :::: name &

# now we have sam files, we can convert it to bam and sort & index them
```
let's check what changed by each step, use N01_AM_Naive.sam as example
here is content without header

```bash
less N01_AM_Naive.sam
SRR7457559.3    4       *       0       0       *       *       0       0       GGGAGGTTCCAGCCAGAGGCTGGAACCTCCCGCATAT   AAAAAEEEAEAEAEEE/EEAEAE/EEEEEEEE<AA/E    YT:Z:UU
SRR7457559.5    0       chr2    27447149        60      75M     *       0       0       CCATTCGGCAGCCTCAAAAACAACAAACCCAACAGGACCTTCTGTTAGACTCTGTATATTATTACTTTTTACAAT      AAAAAEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE     AS:i:0  XN:i:0  XM:i:0  XO:i:0   XG:i:0  NM:i:0  MD:Z:75 YT:Z:UU NH:i:1
```
```bash
mkdir bam
mkdir sorted_bam
mkdir index_bam
nohup parallel -P6 "samtools view -bS sam/{}.sam > bam/{}.bam" :::: name &
nohup parallel -P6 "samtools sort bam/{}.bam -o sorted_bam/{}_sorted.bam" :::: name &
nohup parallel -P6 "samtools index sorted_bam/{}_sorted.bam" :::: name &

```
check bam content, same with sam file
```bash
samtools view N01_AM_Naive.bam | head -n 10

SRR7457559.3	4	*	0	0	*	*	0	0	GGGAGGTTCCAGCCAGAGGCTGGAACCTCCCGCATAT	AAAAAEEEAEAEAEEE/EEAEAE/EEEEEEEE<AA/E	YT:Z:UU
SRR7457559.5	0	chr2	27447149	60	75M	*	0	0	CCATTCGGCAGCCTCAAAAACAACAAACCCAACAGGACCTTCTGTTAGACTCTGTATATTATTACTTTTTACAAT	AAAAAEEEEEEEEEEEEEEEEEEEEE6EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:75	YT:Z:UU	NH:i:1
```

check sorted bam
```bash
samtools view N01_AM_Naive_sorted.bam | head -n 10
# sorted by the start position of the chromsomes
SRR7457559.5779826	272	chr1	3015214	1	75M	*	0	0	TATGGGATGGATCCCTGCATATGGCAATCACTAGATGGTCCATCCTTTTGTCACAGCTCCAAATTTTGTCTCTGT	EAEEEAE//A///</EEEEEEEAAEEEEEA/EAEAEEAAEEEEEEEEEEEEAEEAE<EEAAEEEEEEEEEAAAA/	AS:i:0	ZS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:75	YT:Z:UU	NH:i:5
SRR7457559.2409229	256	chr1	3016283	1	74M	*	0	0	TCGAGGCTTTTCCCTACTTTCTCCTCTGTAAGTTTCAGTGTCTCTGGTTTTATGTGGAGTTCCTTAATCCACTT	6AAAAEAEAEEEEEAEEEEAEEEE/EEEEAAEEEEEEEEEEEEEEAEEEEEEEEE/EEEEEEEEEEEEAEEEEE	AS:i:0	ZS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:74	YT:Z:UU	NH:i:5
```

for indexed bam file, we will get .bai file, which is binary format that we can not directly investigate them.
```bash
#we can use
samtools idxstats N01_AM_Naive_sorted.bam
# to check the alignment statistics
### part of the content 
chr1	195471971	1487358	0
chr10	130694993	1565220	0
chr11	122082543	1848874	0
chr12	120129022	787946	0
chr13	120421639	823018	0
chr14	124902244	834054	0
chr15	104043685	926670	0
chr16	98207768	630159	0
chr17	94987271	1165790	0
chr18	90702639	466556	0
chr19	61431566	1114124	0
chr1_GL456210_random	169725	550	0
chr1_GL456211_random	241735	818	0
###
```
the second column is the length of reference sequence(bp) while the third column is the number of mapped reads to this reference sequence, the fourth column is the number of unmapped reads, which is 0 in our file, indicates there are no unmapped reads in this bam file.

### 4.Quantify gene expression
before we do the quantification, let's check our annotation file

```bash
head Mus_musculus.GRCm38.99.gtf
#!genome-build GRCm38.p6
#!genome-version GRCm38
#!genome-date 2012-01
#!genome-build-accession NCBI:GCA_000001635.8
#!genebuild-last-updated 2019-09
1	havana	gene	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693"; gene_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC";
1	havana	transcript	3073253	3074322	.	+	.	gene_id "ENSMUSG00000102693"; gene_version "1"; transcript_id "ENSMUST00000193812"; transcript_version "1"; gene_name "4933401J01Rik"; gene_source "havana"; gene_biotype "TEC"; transcript_name "4933401J01Rik-201"; transcript_source "havana"; transcript_biotype "TEC"; tag "basic"; transcript_support_level "NA";
```
it records the information of mm10 annotation, havana is the project name, we can ignore for now, second column is feature type, 3rd and 4th column is the start and end position of this feature, . + . mean score, strand and frame, which score and frame is not applicable here. the last column contains gene information which could be useful in the further analyze.

let's do feature count based on bam file, align reads number to each gene id.
this article have the detailed information of featurecounts algorithm: https://academic.oup.com/bioinformatics/article/30/7/923/232889

```bash
nohup featureCounts -T 8 -t exon -g gene_id -a Mus_musculus.GRCm38.99.gtf -o counts.txt sorted_bam/*.bam &
```

let's check out counts.txt
it's a matrix which contains reads of each sample in each geneid (we specify the type of exon, so this gene ids are exons)

```bash
Geneid  Chr     Start   End     Strand  Length  sorted_bam/N01_AM_Naive_sorted.bam      sorted_bam/N02_AM_Naive_sorted.bam      sorted_bam/N03_AM_Naive_sorted.bam      sorted_bam/N04_AM_Naive_sorted.bam
        sorted_bam/R01_AM_Allo2h_sorted.bam     sorted_bam/R02_AM_Allo2h_sorted.bam     sorted_bam/R03_AM_Allo2h_sorted.bam     sorted_bam/R04_AM_Allo2h_sorted.bam     sorted_bam/R05_AM_Allo24h_sorted.bam
        sorted_bam/R06_AM_Allo24h_sorted.bam    sorted_bam/R07_AM_Allo24h_sorted.bam    sorted_bam/R08_AM_Allo24h_sorted.bam
ENSMUSG00000102693      1       3073253 3074322 +       1070    0       0       0       0       0       0       0       0       0       0       0       0
ENSMUSG00000064842      1       3102016 3102125 +       110     0       0       0       0       0       0       0       0       0       0       0       0
```

### 6.further analyze 
now we are ready for using R to conduct further data visualization















