# DRIP-seq learning and analyze notes
### Yifan
## mechanism
DRIP stands for DNA-RNA immunoprecipitation, by directly obtain the sequence of DNA, researchers can further investigate R-loop mapping in the genome.


a monoclonal antibody, normally S9.6 antibody is usd to immunoprecipitate the R-loops, as this antibody could recognizes RNA-DNA hybrids,then non-specific DNA is washed away, leaving RNA-DNA enriched samples. By conducting amplification and standard library preparation, we can using sequencing to capature R-loops information in genome-wide. 


The sequencing data would go through the standard quality-check, alignment process, lastly conduct peak calling to identify regions of significant R-loop formation.

for the wet-lab protocol, this paper provide very helpful information:
Sanz LA, Castillo-Guzman D, Chédin F. Mapping R-Loops and RNA:DNA Hybrids with S9.6-Based Immunoprecipitation Methods. J Vis Exp. 2021 Aug 24;(174):10.3791/62455. doi: 10.3791/62455. PMID: 34515688; PMCID: PMC9676068.
https://pmc.ncbi.nlm.nih.gov/articles/PMC9676068/#S2

## data analysis
for data analyzing, we may link with ChIP-seq concepts, as the inner logic is using cross-linking to capature an interaction, then ampliy, sequencing. 
the step would start from quality control, alignment, peakcalling, and annotation and enrichment.
Let's do together in the next section!


