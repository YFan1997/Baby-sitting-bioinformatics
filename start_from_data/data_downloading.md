# there are many open source for data capturing, here I am going to introduce some popular source.

## NCBI would be one of the most popular website
## directly download by clicking the link is always a choice
## for multiple files downloading,let's use a quick and easy way 
    for example, we need to download several files named:
      SRR111111, SRR12345, SRR123478 .....
    let's open terminal or linux operating system

```bash
## define name and manually writing is a choice, carefully about space, it also could be generated by clicking Accession List button, let's name it as SRA_list, for manually input, use "nano SRA_list", then we can paste or type the asscession index. (use cat SRA_list to inspect)

## use "parallel" for quickly downloading, we will directly have fastq data

parallel 'prefetch {} && fasterq-dump {} -O /path/to/data' :::: SRA_list

## rename fastq data for further analyze
## for example
mv example_1.fastq my_sample_1.fastq

## if many sample, let's use a while loop in terminal, create a txt file by nano which each line is the sample name,

while read -r SRA_list Sample_name; do
  mv "${SRA_list}"*.fastq "${Sample_name}"*.fastq
done < input_file

## in the input file:
SRA123 SampleA
SRA456 SampleB

# what we originally have
SRA123_1.fastq
SRA123_2.fastq
SRA456_1.fastq
SRA456_2.fastq

# after run
SampleA_1.fastq
SampleA_2.fastq
SampleB_1.fastq
SampleB_2.fastq
```
