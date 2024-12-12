# RNA-editing principles and analyze
RNA editing is a post-transcriptional modification of RNA, the transcripts will be edited and make a different transcript, further affecting gene expression. Some common type in RNA editing including, **A-to-I** (adenosine-to-inosine) and **C-to-U** (cytosine-to-uracil) etc. This edition action is mainly conducted by RNA deaminases, such as Apobec1 and ADAR2. RNA-seq will be utilize to capture the sequencing, and variant calling can be applied for processed data.

The usage of RNA-editing is widespreading, here I am going to introduce the one with RNA Binding Protein (RBP), imagine link the RNA-editing enzyme with the target protein, when the protein binding with specific region, the RNA-editing enzyme functions, changes captured, it's acutally indicationg the position of the RBP region.

There are several papers pointed out this method and got promising results:
>McMahon, A. C., Rahman, R., Jin, H., Shen, J. L., Fieldsend, A., Luo, W., & Rosbash, M. (2016). TRIBE: hijacking an RNA-editing enzyme to identify cell-specific targets of RNA-binding proteins. Cell, 165(3), 742-753.

>Ruan, X., Hu, K. & Zhang, X. PIE-seq: identifying RNA-binding protein targets by dual RNA-deaminase editing and sequencing. Nat Commun 14, 3275 (2023). https://doi.org/10.1038/s41467-023-39054-8

>Brannan, K.W., Chaim, I.A., Marina, R.J. et al. Robust single-cell discovery of RNA targets of RNA-binding proteins and ribosomes. Nat Methods 18, 507–519 (2021). https://doi.org/10.1038/s41592-021-01128-0

when compare experiment genome with the reference genome, we can find the difference between the treatment and the reference, just like the SNP. I am going to use PIE-seq as the main reference data and method, which apply two enzymes together.

## preprocessing 
define your path and try the code below
```bash
## downloading the raw data from the GSE155844, here I only select several 
parallel 'prefetch {} && fasterq-dump {} -o raw.data' :::: SRR_list 

#rename input_file by mv

#nano input_file
SRR14612788 empty_1
SRR14612789 empty_2
SRR14612790 APAD_48_1
SRR14612791 APAD_48_2
SRR14612792 PUM2_48_1
SRR14612793 PUM2_48_2

while read -r SRR_list Sample_name; do
  for fastq_file in raw.data/${SRR_list}*.fastq; do

    # Extract the suffix (_1 or _2, etc.) from the original file

    suffix=$(echo "$fastq_file" | sed "s/.*${SRR_list}\(.*\)\.fastq/\1/")
    
    # Rename the file to match the new Sample_name
    mv "$fastq_file" "raw.data/${Sample_name}${suffix}.fastq"
  done
done < input_file

## in the next code, the {} stands for the filename, I recommend using one file for first try, if works, apply parallel, but be careful of the computational capability.

## in their method supplemental documentation, it list the parameter used for STAR alignment, so just apply here
STAR --twopassMode Basic --runThreadN 16 \
--genomeDir /mnt/pv_compute/yifan/ref_genome/index/hg38_star \
--genomeLoad NoSharedMemory \
--outFileNamePrefix $path/STAR_refine_parameter/{}_ \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--outSAMstrandField intronMotif \
--outSAMattributes All \
--outFilterType BySJout \
--outFilterMultimapNmax 1 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--readFilesIn {}_1.fastq {}_2.fastq
```

## Remove duplication
is very important in RNA-editing analyze to remove the contanmination caused by PCR, and make the result more clear.
the tool used is picard, which need java to implement it (most of computer already have java)
download from 
https://broadinstitute.github.io/picard/

then $picard="java -jar $path/picard.jar"  (carefull about the version)

the tool of variant calling is JACUSA2, https://github.com/dieterich-lab/JACUSA2, using by the same flow: download and java.

jacusa2="java -jar $path/JACUSA_v2.0.4.jar"
```bash

## remove duplicate
$picard MarkDuplicates \
I= {}_Aligned.sortedByCoord.out.bam \
O= $path/{}_dup.bam \
M= {}.rmdup.log \
REMOVE_DUPLICATES= true \
ASSUME_SORTED= true \
VALIDATION_STRINGENCY=LENIENT

## index bam
samtools index {}_dup.bam

## vairant calling example, need to apply for several groups
## Empty vs APAD 
## Empty vs PIE
## APAD vs PIE

$jacusa2 call-2 -s -c 5 -P RF-FIRSTSTRAND -p 20 -W 1000000 -F 1024 \
-filterNM_1 7 -filterNM_2 7 -T 1 -a D,Y -r Empty.vs.ADAR.all.editing.sites \
$path/empty_1_DupRm.bam,$path/empty_2_DupRm.bam \
$path/APAD_48_1_DupRm.bam,$path/APAD_48_2_DupRm.bam

## intersect with common SNP 
bedtools intersect -v -a Empty.vs.PUM.all.editing.sites -b dbsnp_138.hg38.vcf.gz > Empty.vs.PUM.DP5.wo.snp

## REMEMBER TO RE-ADD HEARDER for further analyze, seperate by tab
#contig start   end name    score   strand  bases11 bases12 bases21 bases22 info    filter  ref
```
The Markdown file is how to use JACUSA2 helper to further analyze the result, combine with GO and KEGG. apply to your processed result!

