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
nohup parallel "hisat2 -p 8 -x mm10_index -U {}.fastq -S sam/{}.sam" :::: name &












