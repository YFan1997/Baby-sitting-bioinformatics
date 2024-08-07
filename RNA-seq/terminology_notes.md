## concepts and type of analysis

### 1. R packages

we actually mainly use DESeq2, but other packages are also very useful for deeper understanding and management of our result.
In DESeq2, there are many embeded algorithms which users can directly use, but, I really want to digging in the meaning of these parameters, we will find these would help a lot in understanding plots and make our own plots based on specific needs.

### 2.Understanding some concepts
> log2FC


this parameter stands for log2 fold change. firstly, fold change measures how much quantity changes between two conditions.let's say condition A in gene x's expression level(read count) is a, and condition B in gene x' expression level(read count) is b, then the fold change will be b/a, the log2FC will be log2(b/a). Why using log2, as it could make data symmetric, as if b = a, log2(b/a) = 0, etc.  this concept will help us understanding plotMA

> pvalue and padj

in the result generated from DESeq package, we can extract the result easily by results() function. the pvalue infers the statistical significance of the observed difference in expression for a gene between conditions being compared. usually used Wald test; the padj is adjusted pvalue, when we testing thousands of genes, some may just appear siginificant by chance, in DESeq2, it uses Benjamini-Hochberg procedure to adjust p-values, which ranks p-values and adjust to control false discorvery rate. so in some circumstances, padj is more reliable, and we may set threshold for both pvalue and padj.



> data transformation: vst, rlog, ntd

this step could help visualizating and clustering.  when we have a RNA-seq data, the read counts exhit a strong dependence of the variance, this may due to the data generation process itself，as RNA-seq counts are normally modeled by Poisson-distributed, which the variance is equal to the mean, for lowly expressed genes, the variance is low because the mean count is low, so as highlt expressed genes, the variance increases as the mean count increase. VST and rlog stand for variance stabilizing transformation and regularized logarithm. Both of them tend to remove the dependence of variance on the mean.


> normTransform in DESeq


this is a built in function performs normalization by transforming the raw counts to the normalized counts which use size factor estimated by DESeq2, without applying variance-stablizing transformation, which could be useful for creating heatmaps and PCA plot.


When using meansdplot to get the plot, the system will caculate the mean counts of each gene, and rank them by descending, this formed x-axis, each gene' standard deviation(sd) formed y-axis, also the density will performed, which indicates how many gene share the similar mean and sd.
after conduct vst, what should we expect? Even though the real data can be very chaos, if we see a relatively flat trendline, which means VST has done its job to migate the high variance in high mean expression values.


>sample-to-sample distance

this parameter is to calculate the elucidean distance between samples, after obtain vst result, transpose the matrix, which makes sample as row, while the gene as the column, then conduct calculation between samples to check the similarity between samples, and small distance infers high similarity.


> PCA

it stands for principle component analysis, which is based on mathematical calculation based on covariance matrix, thankfully, the embeded programming do all these calculation for us, but it's always good and worth to understand what the programming is doing. In summary, it peaks the most siginificant PC to plot the graph, in most cases, the PC1, PC2 capture most of the variantion in the data. When we get the PCA plot, the similar sample tend to clustring together. Very easy for us to interpret the plot.

## 3.Understanding analysis
> GO analysis


GO stands for gene ontology, which is powerful for categorize and interpret the roles of genes. There are several components of gene function : BP (biological processes), CC (cellular components), MF (molecular functions). It's an enrichment analysis.



> GSEA analysis


after conducting gseGO, the enrichment score (ES) will be calcuted, which is reflect the degree of a gene being overrepresent, for simply understanding, positive score means upregulated, which can be called activated gene, wthile the negetive score means downregulated, which can be called suppressed gene. NES stands for normalized enrichment score, which conduct a normalization through gene sets for comparing. 

> KEGG analysis


KEGG stands for kyoto encyclopedia of genes and genomes, which contains several database to help understaning pathway, protein information, when we doing KEGG analysis, we actually doing KEGG gene enrichment analysis, to identify active biological process based on these database annotation. 

> Pathview


This method is used to visualize gene or pathway from KEGG pathways. we will take a look at it later



