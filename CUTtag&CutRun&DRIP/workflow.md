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
## unzip by gunzip mm10.fa.gz
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

## remove duplicate
mkdir -p result/spike_rmdup
mkdir -p result/rmdup

#picard="java -jar path/to/picard.jar"

parallel "$picard MarkDuplicates \
I= spikebam/{}_sorted.bam \
O= result/spike_rmdup/{}_rmDup.bam \
M= result/spike_rmdup/{}.rmdup.log \
REMOVE_DUPLICATES= true \
ASSUME_SORTED= true \
VALIDATION_STRINGENCY=LENIENT" :::: samplenames

parallel "$picard MarkDuplicates \
I= bamout/{}_sorted.bam \
O= result/rmdup/{}_rmDup.bam \
M= result/rmdup/{}.rmdup.log \
REMOVE_DUPLICATES= true \
ASSUME_SORTED= true \
VALIDATION_STRINGENCY=LENIENT" :::: samplenames

## normalization
## in this step, we will use spike in reads to calculated the scaling factor, please note here, there is no unique way to calculate scaling factor, here I am going to use
 mapped reads of dm6 /mapped reads of mm10

## get total reads number

parallel "echo {}: \$(samtools view -c result/rmdup/{}_rmDup.bam)" :::: samplenames
parallel "echo {}: \$(samtools view -c result/spike_rmdup/{}_rmDup.bam)" :::: samplenames


# here I got total read
## for mm10:
#KO_rep1: 13797090
#KO_rep2: 15004470
#KO_rep3: 16468037
#WT_rep2: 19908553
#WT_rep1: 22842614
#input_WT: 25471756
#WT_rep3: 29826345
#input_KO: 29572179

## for dm6

#KO_rep1: 23132552
#KO_rep2: 26471484
#KO_rep3: 29499876
#WT_rep2: 29251105
#WT_rep1: 31165283
#input_WT: 31080133
#input_KO: 36044150
#WT_rep3: 42200140

## get mapped read
parallel "echo {}: \$(samtools view -c -F 4 result/rmdup/{}_rmDup.bam)" :::: samplenames

parallel "echo {}: \$(samtools view -c -F 4 result/spike_rmdup/{}_rmDup.bam)" :::: samplenames
## I obtained
## mm10
#KO_rep1: 13183342
#KO_rep2: 14343264
#KO_rep3: 15610472
#WT_rep2: 19108384
#WT_rep1: 21848101
#input_WT: 25173676
#WT_rep3: 28497478
#input_KO: 29233638

## dm6
#KO_rep1: 152398
#KO_rep2: 163492
#KO_rep3: 183109
#WT_rep2: 226949
#WT_rep1: 291376
#input_WT: 114706
#input_KO: 128248
#WT_rep3: 297802

## we can use excel sheet to drag the scaling factor out
## scaling factor

#input_KO:	0.004387001
#input_WT:	0.004556585
#KO_rep1:	0.011559891
#KO_rep2:	0.011398521
#KO_rep3:	0.011729882
#WT_rep1:	0.013336445
#WT_rep2:	0.011876933
#WT_rep3:	0.010450118

## convert bam to bed
mkdir -p result/bed

# get mapped bam read to bed file, note here the read is single-read, so use bamtobed -i -; if is paried, used -i -bedpe

nohup parallel "samtools view -b -F 4 bamout/{}_sorted.bam | bedtools bamtobed -i - > result/bed/{}.bed" :::: samplenames &

# while it running, perpare scaling.factor file by nano, and paste the result from previous step, replace semicolon with tab

# Read scaling factors and apply them
mkdir -p result/normalized_bed

while read sample factor; do
    awk -v scale=$factor 'BEGIN {OFS="\t"} {if (NF >= 5) $5 = $5 * scale; print}' \
        result/bed/${sample}.bed > result/normalized_bed/${sample}_normalized.bed
done < scaling.factor

## let's have a check
head input_KO_normalized.bed
chr1	3000090	3000130	SRR14320692.32728410	0	+
chr1	3000117	3000191	SRR14320692.17471565	0.171093	+
chr1	3000299	3000374	SRR14320692.24556080	0.179867	+
chr1	3000324	3000399	SRR14320692.17107622	0.241285	-
chr1	3000485	3000558	SRR14320692.30422131	0.04387	+
chr1	3000622	3000696	SRR14320692.4006276	0.065805	-
chr1	3000870	3000945	SRR14320692.11132755	0.26322	+
chr1	3000927	3001001	SRR14320692.31862647	0.114062	-
chr1	3000981	3001049	SRR14320692.3858590	0.26322	+
chr1	3001022	3001097	SRR14320692.19683364	0.26322	+

head input_KO.bed 
chr1	3000090	3000130	SRR14320692.32728410	0	+
chr1	3000117	3000191	SRR14320692.17471565	39	+
chr1	3000299	3000374	SRR14320692.24556080	41	+
chr1	3000324	3000399	SRR14320692.17107622	55	-
chr1	3000485	3000558	SRR14320692.30422131	10	+
chr1	3000622	3000696	SRR14320692.4006276	15	-
chr1	3000870	3000945	SRR14320692.11132755	60	+
chr1	3000927	3001001	SRR14320692.31862647	26	-
chr1	3000981	3001049	SRR14320692.3858590	60	+
chr1	3001022	3001097	SRR14320692.19683364	60	+

# which 0.171093 = 39 * 0.004387001

## let's do peakcalling
## here we using SEACR for narrow peak, which require bedgraph format
##https://github.com/FredHutch/SEACR
## the other tools can also be used, such as MACS2

mkdir -p result/peakcalling
## convert bed to bedgraph
## mm10 chromsize:
wget https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes
mkdir -p result/bedgraph

## make sure the fragement is smaller than 1000
parallel "awk '\$3-\$2 < 1000 {print \$0}' result/normalized_bed/{}_normalized.bed > result/bedgraph/{}.clean.bed" :::: samplenames

## peak chr, start end signal information and sorted
parallel "cut -f 1,2,3,5 result/bedgraph/{}.clean.bed | sort -k1,1 -k2,2n -k3,3n > result/bedgraph/{}.fragments.bed" :::: samplenames

## convert bed to bedgraph based on mm10 chromosize
parallel "bedtools genomecov -bg -i result/bedgraph/{}.fragments.bed -g mm10.chrom.sizes > result/bedgraph/{}.fragments.bedgraph" :::: samplenames

## apply SEACR
$seacr="path/to/SEACR_1.3.sh"

bash $seacr result/bedgraph/KO_rep1.fragments.bedgraph \
result/bedgraph/input_KO.fragments.bedgraph \
non stringent result/peakcalling/KO_rep1_vs_control.peaks

## MACS2 code ( the reason I choose MACS2 as SEACR capture very limited peaks, as no software is perfect, we can alway consider about other software)
## to more directly see the peak, use IGV to load bed files
##https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html

macs2 callpeak -t KO_rep1.clean.bed -c input_KO.clean.bed -f BED -g mm -n KO_rep1_vs_control

## Seacr works well in wild type
bash $seacr result/bedgraph/WT_rep1.fragments.bedgraph \
result/bedgraph/input_WT.fragments.bedgraph \
non stringent result/peakcalling/WT_rep1_vs_control.peaks

## As I am using these two software, even though they both be implemented with good algorithms, the SEACR requires additional data handling
## and get different result from macs2, depends on situation, choose carefully.

## I am sharing my shiny.app for analyze bed file, which could directly upload bed file and get annotation (so far is only for hg38, so not suitable for this data)
https://babysittingbioinfo.shinyapps.io/shiny_r_project/
```
