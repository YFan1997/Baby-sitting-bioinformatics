# CUT&RUN, CUT&Tag, and DRIP-seq

In this section, rather than focusing on a single method, I will discuss these three techniques together to highlight their differences and similarities. I will also introduce how to analyze their raw data and interpret the results. Since all of these methods are conceptually based on ChIP-seq, if you are familiar with one, you will find it straightforward to work with the others. Otherwise, once you learn one method, you can easily apply the same principles to all of them.

## Mechanism

**CUT&RUN (Cleavage Under Targets and Release Using Nuclease)**  
In CUT&RUN, an antibody specific to a protein of interest is used in combination with a DNA-cleaving enzyme (often micrococcal nuclease). The enzyme binds to the antibody-protein-DNA complex and cuts the DNA near the labeled sites, releasing short DNA fragments.

**CUT&Tag (Cleavage Under Targets and Tagmentation)**  
CUT&Tag employs a modified Tn5 transposase that can insert sequencing adapters (tagmentation) directly into the target DNA. This approach streamlines the library preparation process by simultaneously fragmenting and tagging the DNA, thus saving time.

**DRIP (DNA-RNA Immunoprecipitation)**  
DRIP targets R-loop structures using the S9.6 antibody, which specifically recognizes RNA-DNA hybrids. After binding, these complexes are immunoprecipitated, enriching for DNA regions that contain R-loops.

### Commonality and Peak Calling

All three methods rely on antibodies to enrich specific genomic regions, ultimately generating DNA sequencing data. The analysis involves identifying "peaks" where these antibodies have enriched certain genome intervals. Peak calling, a concept widely used in ChIP-seq, is also applied here to determine regions of significant enrichment.

### Additional Resources

For additional information, the following articles provide valuable information, especially the CUTTag_tutorial, provide a detailed example of how to analyze CUTTag data:

- Sanz LA, Castillo-Guzman D, Chédin F. *Mapping R-Loops and RNA:DNA Hybrids with S9.6-Based Immunoprecipitation Methods.* J Vis Exp. 2021 Aug 24;(174):10.3791/62455.  
  [https://pmc.ncbi.nlm.nih.gov/articles/PMC9676068/](https://pmc.ncbi.nlm.nih.gov/articles/PMC9676068/)

- https://www.cellsignal.com/applications/cut-tag-overview
- https://yezhengstat.github.io/CUTTag_tutorial/

## Data Analysis


For data analysis, you can follow ChIP-seq concepts since the underlying logic is similar. Once you’re familiar with analyzing one method, you can apply the same principles—quality control, alignment, and peak calling—to the others.

In this tutorial, I am going to using DRiP seq data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173258 which contain two input data for Knockout model, and Wide type, which can serve as control; 6 treatment data, among three are KO model, three are WT model.

```bash
##rename_map
SRR14320691 input_WT
SRR14320692 input_KO
SRR14320693 WT_rep1
SRR14320694 WT_rep2
SRR14320695 WT_rep3
SRR14320696 KO_rep1
SRR14320697 KO_rep2
SRR14320698 KO_rep3

parallel "prefetch {} && fasterq-dump {}" :::: srrlist

while read old new; do
    mv ${old}.fastq ${new}.fastq
done < rename_map

for file in *.fastq; do
    basename "$file" .fastq
done > samplenames

## I found this step may create one summary line in header, for clean, copy the samplenames and paste

## Quality control
parallel "fastqc {}.fastq -o fastqc_output/" :::: samplenames


## trim
# trimmomatic resource
#http://www.usadellab.org/cms/?page=trimmomatic
trimmomatic="java -jar /mnt/pv_compute/yifan/practice/temp/Trimmomatic-0.39/trimmomatic-0.39.jar"
# perform quality trimming using sliding window of 5 nucelotids, and if average quality in the window drops below Q15 trimming occurs. and discar the reads shorter than 15 nucleotides 
parallel "$trimmomatic SE -phred33 {}.fastq output.trimmed_{}.fastq SLIDINGWINDOW:5:15 MINLEN:15" :::: samplenames

## after trimming. we can use fastqc check again.

## alignment
## Using bwa for alignment
## first build index for mm10. mm10.fa can be downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
## wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
## for dm6: wget https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz
## unzip by gunzip mm10.fa.gz etc.
## the dm6 is spikein genome, will be used for further normalization.


mkdir samout
mkdir spikesam
mkdir bamout
mkdir spikebam

## here I am using nohup for hang the command in the server, so you won't worry about if your computer closed
## you can check status by "jobs"
## the following steps take times

bwa index mm10.fa
bwa index dm6.fa
nohup parallel "bwa mem -t 8 /mnt/Data1/HwangLab/Yifan/ref_genome/mm10.fa trimmed_data/output.trimmed_{}.fastq > samout/{}.sam" :::: samplenames &

nohup parallel "bwa mem -t 8 /mnt/Data1/HwangLab/Yifan/ref_genome/dm6.fa trimmed_data/output.trimmed_{}.fastq > spikesam/{}.sam" :::: samplenames &

## convert sam to bam

parallel "samtools view -bS samout/{}.sam > bamout/{}.bam" :::: samplenames
parallel "samtools view -bS spikesam/{}.sam > spikebam/{}.bam" :::: samplenames

## sort and index bam file
parallel "samtools sort bamout/{}.bam -o bamout/{}_sorted.bam" :::: samplenames
parallel "samtools index bamout/{}_sorted.bam" :::: samplenames

parallel "samtools sort spikebam/{}.bam -o spikebam/{}_sorted.bam" :::: samplenames
parallel "samtools index spikebam/{}_sorted.bam" :::: samplenames
```
