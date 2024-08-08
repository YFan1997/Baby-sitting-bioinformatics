# DNA methylation
## Bisulfite treatment 
In the experiment, we can use bisulfite treatment to detect the location of DNA methylation. The mechanism is that bisulfite treatment does not affect methylated cytosine (C-methy), but it converts unmethylated cytosine (C) to thymine (T). After sequencing, we can determine the methylation status by checking whether a base remains as cytosine (indicating methylation) or has been converted to thymine (indicating lack of methylation)

## data description 
### illumina DNA methylation arrays
after bisulfite treatment, the bisulfite-converted DNA will be fragmented and amplified using PCR. These amplified DNA will be hybridized to the array, which contains thousands of probes specific to CpG sites, and then conducted single-base extension for adding a fluoresently labled single nucleotide to the hybridized probe, normally, a green fluorescent tag would be used for the methylated C, red for unmethylated U. Then the scanner would detect the fluorescent signals at each probe site, these signals can be read as raw intensity data.

> methylation beta values
this term refer the methylation level at each CpG site, range from 0 to 1, which can be calculated by Beta = M/(M+ U + α)  which M is the signal intensity of the methylated data, U is the signal intensity of the unmethylated probe, α is a constant to prevent division by zero. In R, we can simply use functions to get this value.

### Whole genome Bisulfite Sequencing
this is another way for obtaining DNA methylation data, during the sequenceing, the DNA fragments are read in a massively parallel manner and generating large amounts of data including original sequence and bisulfite-converted bases. Normally, we use specific alignment tools such as Bismark to align the result in reference genome. After alignment, the methylation status of each C will be determined, and the proportion of reads C/T will be calculated.

after obtaining these based information, we can use R to analyze real-world data.




