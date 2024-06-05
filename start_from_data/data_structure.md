# This page would talks about fasta,fastq,sam and bam in baby-sitting way
## fasta
this format is a simplified sequence storage file, end with .fasta or fa, normally, it's constructed by
```bash
>sample_name1|sample_id1
ATGCTAGCATCGAGCATCGACATCGACTAGCATCACTAGC....\
>sample_name2|sample_id2
ATGCTAGCATCGAGCATCGACATCGACTAGCATCACTAGC....
.....

# we can use grep in linux to obtain specific id information
grep ">" file1.fasta

# some useful linux command, more to be introduced in code tips folder
less file1.fasta # check the content of file, useful for checking big file as it shows in pages
grep ">" file1.fasta | wc -l #check how many header information in fasta file
```

## fastq
this format extensd FASTA with quality of sequencing and more information in header, sometimes we would have data
sample1_1.fastq.gz  sample1_2.fastq.gz, this indicates that the sample is sequenced in pairend method
```bash
@A00695:90:HNLJVDRX2:1:2101:7310:1000 1:N:0:NAAGGCGA+NCGATCTA
NGAAGCAAAAGCAGAAACCCCCGATAAACCCATTGCATTTTGTGAAACTTATTCCCTATAACAAGAGTAGCACGGGAAAGACTGGCCCCCATGAGTCAATTACCTCCCCCTGGGTCCCTCCCACAACATGTGGAAATTCTGGGACATACA
+
#FFF,FFF,,FFFFFFFF,F,F:FFFFF:,FFF:,:FFFF,FF,FFF:FFFFF:,:,FF,FFFF:FF:F,FFF,FFFFFFFF:F,FFFF,F,,F:,FFF:,F,,FFF:,::F,FFF:FF:FFFFFFFFF,,::FFF,::F,F,F,F::,F

# the line start with "#" is the quality metrix for each base, there is a table in this website:https://help.basespace.illumina.com/files-used-by-basespace/quality-scores

```

## sam
there are two main sections in sam format of data, one is header, another one is alignment section. In header
@SQ stores the information about the reference sequences, LN stores the sequence length. SN stands for sequence name.
@PG contains metadata about the programs used to generate SAM or BAM files.
the line without start with @ is the start of alignment section

```bash
# check the file
samtools view sample.sam
```
in the alignment section, it contains query name, flag, reference name, position on the reference sequence, mapping quality, CIGAR string, Reference name and position, length, original read sequence and orginal read base quality.
Among these information, Flag and Mapping Quality can be used to infer important sequencing quality
for example if we have a flag equal to 83, we can convert it to binary by python
```python
print(bin(83))
```
which will give us 0b1010011, we can ignore the b, consider its as 9 digit string, 1 for true, 0 for false,
stare from right as position 1, and then in linux,
```bash
samtools flag
```
to check each position's meaning, for example, the position 1 is 1, so it's a pair end sequence. etc...
we can directly check by
```bash
samtools flags 83
```
## bam
we can use samtools to conver sam file to the bam file, bam file can be considered as a binary version of sam file.
```bash
samtools view -b sample.sam > sample.bam # for include header, using -h
# other common linux command including
# samtools index
# samtools sort
```
now we have covered the basic data structure in bioinformatics, you are ready for further practical content.


