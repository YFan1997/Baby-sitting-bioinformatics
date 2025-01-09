# Hi-C

## Mechanism
Hi-C is an advanced technique based on the 3C (Chromosome Conformation Capture) method. While the previous sections covered DNA-RNA and DNA-protein interactions, this chapter focuses on another important aspect of molecular interactions: chromatin interactions.

Chromatin is composed of DNA and histone proteins, with nucleosomes serving as its fundamental repeating units. Formaldehyde, a widely used crosslinking reagent, is applied in this context to fix DNA and proteins in their native spatial arrangements within the cell. This crosslinking preserves interactions between DNA and proteins that are in close physical proximity. Such spatial arrangements are critical in biological processes, such as the formation of topologically associating domains (TADs), which play a key role in genome organization and gene regulation.

Although both 3C and ChIP (Chromatin Immunoprecipitation) techniques rely on DNA-protein interactions, their methodologies differ significantly. Unlike ChIP, which requires specific antibodies to target particular proteins, 3C does not depend on antibodies. Instead, it leverages the natural interactions between DNA and proteins to capture the spatial organization of chromatin. For example, mechanisms like transcription factories bring distant DNA regions into close proximity for coordinated gene expression.

Hi-C extends the 3C technique to a genome-wide scale. Following restriction enzyme digestion, biotin is used to label the ends of DNA fragments. These fragments are then ligated together under dilute conditions, capturing interactions between DNA regions that were spatially close. Unligated fragments are removed, and streptavidin-coated beads are employed to isolate the biotin-labeled ligation products. The purified DNA is then ready for sequencing, enabling comprehensive analysis of chromatin interactions across the genome.

Reference:Jon-Matthew Belton, Rachel Patton McCord, Johan Harmen Gibcus, Natalia Naumova, Ye Zhan, Job Dekker,Hi–C: A comprehensive technique to capture the conformation of genomes, Methods, Volume 58, Issue 3, 2012,Pages 268-276, ISSN 1046-2023,https://doi.org/10.1016/j.ymeth.2012.05.001.

helpful video illustration: https://www.youtube.com/watch?v=F6TRAqKkgoo 
https://www.youtube.com/watch?v=0787ciSyrT8

## Data
the valid reads in HI-C is from different restriction fragments and facing toward a restriction site. 










