
# Downstream Analysis of RNA-Seq Results in R: GSEA, PPI Networks, and Biological Interpretation

## Workshop Overview
After identifying differentially expressed genes as we did in the [previous workshop]((https://github.com/merlinis12/RNA-Seq-Data-Analysis-in-R/wiki)), it’s critical to interpret these findings in a biological context. 
This workshop explores downstream analyses you can perform with differential gene expression analysis results in R.
You’ll learn how to:
- Conduct [Gene Set Enrichment Analysis (GSEA)](#2-gene-set-enrichment-analysis-gsea) for pathway-level insights.
- Build [Protein-Protein Interaction (PPI) networks](#3-protein-protein-interaction-ppi-networks) to identify functional connections between genes.
- Explore drug repurposing opportunities based on RNA-seq data.


---

## 1. Introduction to Downstream Analysis
Differential gene expression (DGE) analysis is an essential step in RNAseq downstream analysis. The goal is to identify differentially expressed genes (DEGs) between two conditions.
For example, in the [previous workshop](https://github.com/merlinis12/RNA-Seq-Data-Analysis-in-R) we studied the difference in gene expression between airway smooth muscle cells of healthy individuals and airway smooth muscle cells in individuals treated with dexamethasone, a strong synthetic glucocorticoid..

But our DGE analysis returned a long list of differentially expressed genes.

Now, you may have many questions like:
- How do we even start interpreting this?
- Is there a way to summarise this long list of genes and interpret hundreds of DGEs at once?

No worries! The good new is that there is a common approach called **pathway enrichment analysis (PEA)** that allows you to summarize
the long gene list to a shorter and more easily interpretable list of pathways.

>  [!IMPORTANT]
> **Pathway enrichment analysis summarises the long gene list to a shorter and more easily interpretable list of pathways.**
>
> So instead of having a list of more than 2.000 genes, you may get a list of 50 biological pathways. And then, you can check which genes are behind these pathways.



BUT  How does pathway enrichment analysis work?

Pathway enrichment analysis needs 3 ingredients.

First, of course, your gene list of interest for example, a list of differentially expressed genes which you want to summarise.

Second, a list of background genes – for example, all of the genes in the human genome.

Finally, it will take a lists of gene sets. Gene sets are basically groups of related genes. Of course, for the algorithm to know if your list has a lot of genes related to breast cancer, or apoptosis, or cellular respiration, you need to tell it which genes are actually involved in breast cancer, apoptosis, and cellular respiration.

You will find out a bit more about each component later on.

** PUT image of list od overrepresented pathways (https://biostatsquid.com/pathway-enrichment-analysis-explained/)**

PEA essentially compares your gene list to the background list to check if there are certain pathways overrepresented.

**In other words, is our list of differentially expressed genes enriched with genes involved in IL-6 synthesis pathway?**

To answer this question, we can build a contingency table. This will help us determine whether the fraction of genes of interest in the pathway is higher compared to the fraction of genes outside the pathway (so,  background set).

You can have a look at the table below.

We have a column for differentially expressed and a column not differentially expressed genes, and then two rows, for genes that are annotated as being involved in IL-6 production and genes that are not involved in IL6 production.

To simplify things a lot, we will just look at 30 genes. 15 differentially expressed genes were identified and of those, 12 genes were associated with the GO term interleukin-6 production.

We find that 12 out of our 15 differentially expressed genes are involved in interleukin-6 production. We could quite confidently say that our gene list is enriched with genes involved in IL-6 production.

Obviously, we need an objective statistical test to determine what is enriched and what is not.

There are many methods out there, but the one that is commonly used in pathway enrichment analysis is Fisher’s exact test.
Like with most statistical tests, you will obtain a p-value.

If you p-value is really low, you can safely say that your list is overrepresented with genes involved in IL-6 production, in other words, IL-6 production is an important pathway in alcoholic liver disease compared to healthy liver cells.

> [!IMPORTANT]
> Careful! A low p-value does not mean that pathway is upregulated! You would have to check the genes that are actually overrepresented in your list, and see if they are positive or negative regulators of that pathway.


> [!NOTE]
> So this is what pathway enrichment is all about.
>
> You summarise a long list of genes to a shorter list of pathways, with their p-values.

Obviously, it does it not with one, but you with thousands of pathways.

And this bring us to a BIG problem.

The big problem is called multiple testing.

Basically, because we are repeatedly testing a lot of pathways, some pathways will get apparently significant p-values just by chance.

So we might get results that are a bit unexpected… or that just don’t make any sense.

Thank goodness there is a solution for this.

We just need a multiple-testing correction method. The most commonly used method is the Benjamini-Hochberg (BH) correction.

SO IN SUMMARY
>[!NOTE]
>Pathway enrichment analysis takes your gene list of interest and compares it to a list of background genes to check if there are certain pathways that are overrepresented.
>
>So it checks the fraction of your genes annotated to a specific Gene Ontology (GO) term. Then it checks the proportion of genes in the whole genome (your background set) that are annotated to that GO Term.
>
>Then, it gives you a p-value which tells you what is the probability that that pathway is actually overrepresented in your gene list and it wasn’t just coincidence.
---
**Gene sets** = groups of related genes

There are many databases of gene sets out there.

Some of the most widely used are Gene Ontology GO, KEGG or Reactome.

- GO basically focuses on biological processes
- KEGG is more focused on metabolic pathways
- Reactome is a curated database of human molecular pathways.

**Gene Ontology GO**
 [Gene ontology](https://www.geneontology.org/) provides a controlled vocabulary for describing:
 1. Biological Processes (BP ontology): Describes the larger biological pathways or events that a gene product is involved in, such as "cell cycle" or "apoptosis".
 2. Molecular Functions (MF ontology): Represents the specific molecular activity of a gene product, like "protein kinase activity" or "DNA binding".
 3. Cellular Components (CC ontology): Indicates where within a cell the gene product is located, such as "plasma membrane" or "nucleus". 

> [!NOTE]
> Hierarchical structure:
> GO terms are organized in a hierarchical structure where more general terms are parents to more specific child terms. 


> [!NOTE]
> ** It is essential to choose well your background genes to your experiment. **
> There are gene sets for diseases, which gives you groups of genes associated to different diseases, tissues, which gives you groups of genes expressed in specific tissues, transcription factor targets, which tells you which genes are the target of different transcription factors


---

**Gene List**
Results from gene expression analysis often look like this. Some genes are downregulated, some are upregulated, some don’t change. Some changes are significant, and some are not significant at all.
If you just use this list for pathway enrichment analysis, it will not take into account all that information. It will also match genes that are not even differentially expressed.

So you need to first filter your results, by significance and fold change, to keep only differentially expressed genes. The genes of interest.

> [!WARNING]
> Of course, the results can change a lot depending on the cut-offs you set to say what is a DEG and what is not.


Is there a more objective or unbiased way of doing this?

For example, by taking into account the significance and strength of the changes, which is information we have anyway?

Of course there is!

The method is called Gene Set Enrichment Analysis and it a great pathway enrichment analysis method which helps us solve this problem. Find more about it in this other post.

--
OPPURE

Pathway enrichment analysis methods take a list of differentially expressed (DE) genes as input, and identify the sets in which the DE genes are over-represented or under-represented.

So they basically summarise long lists of genes into a shorter list of pathways.

The significance of each pathway is measured by calculating the probability that the observed number of DE genes in a given pathway are simply observed by chance. Lower p-values means that the pathway is actually overrepresented and it was not just by chance.

These approaches are known as Over-Representation Analysis (ORA). + picture

The problem in overrepresentation analysis methods is that we first need to select differentially expressed genes.

> [!CAUTION]
> Overrepresentation Analysis Methods depend a lot in the criteria used to select differentially-expressed genes.


How did they eliminate this dependency on gene selection criteria? Basically, by taking all genes into consideration.

But they do it smartly, because they are not only looking for significant changes in sets of functionally related genes, but also genes with large expression changes.

One of the most popular approaches is Gene set enrichment analysis, or GSEA.
---

## 2. Gene Set Enrichment Analysis (GSEA)

> [!IMPORTANT]
> The steps to perform Gene Set Enrichment Analysis are very similar to overrepresentation analysis methods. The big difference is that the input is not a list of genes, but a ranked list genes.

This basically means that genes are ranked by some score.

A common way of ranking genes is by level of differential expression. The p-values tell us how significant the change is. The log2fold changes tell us the direction and strength of the change, basically if they are upregulated or downregulated.

We can combine both to get a ranked list of genes. This will order the genes both by significance and the direction of change.

> Ranking = Sing(log2FoldChange)*-log10(p-value)


The result is basically the same as with overrepresentation analysis.

We will obtain a p-value, which we need to correct for multiple testing since we are repeatedly testing thousands of gene ontology terms.

This way, we were able to reduce our long list of genes into a more manageable list of biological pathways.

> [!NOTE]
> Another note on rankings… You can just use the log2(fold-change) as a ranking – but it will not take into account genes with a large fold change but non-significant. However, if you have already selected your significant genes, this may be a good option for you.

> So as you can see, gene set enrichment analysis has the advantage that you don’t filter out genes prior to the analysis, and also it takes into account how significant the changes are, and in which direction.

> [!TIP]
> Gene set enrichment analysis takes your ranked gene list of interest and compares it to a list of background genes. By statistically testing the distribution in your list, it determines which pathways are overrepresented.


> [!NOTE]
> there is a new generation of pathway enrichment analysis Topology-based (TB) methods which also take into account dependencies and interactions between genes.
---

Gene Set Enrichment Analysis (GSEA) is a method used to determine whether a predefined set of genes shows statistically significant differences in expression between two biological conditions. GSEA does not rely on a threshold (like fold-change) and instead assesses whether the genes in a set are overrepresented at the extremes of a ranked list of genes.




**Steps in GSEA**

- **Rank the genes** based on their differential expression.
- **Predefine gene sets**: These could be biological pathways, gene ontology (GO) terms, or other gene collections.
- **Calculate enrichment scores** for each gene set.
- **Assess statistical significance** using permutation testing.


### 2.1 Overview
GSEA helps identify pathways and biological processes enriched in a ranked list of genes. We’ll analyze pathways using:
- **Gene Ontology (GO)**: Biological processes, molecular functions, and cellular components.
- **KEGG**: Metabolic and signaling pathways.
- **Reactome**: High-level pathway annotations.

### 2.2 GSEA in R with `fgsea`
`fgsea` package in R provides a fast and efficient method for performing preranked GSEA. It enables precise and rapid calculation of extremely low GSEA p-values across a collection of gene sets. The p-value estimation employs an adaptive multi-level split Monte Carlo approach. For a detailed explanation of the algorithm, refer to [this preprint](https://www.biorxiv.org/content/10.1101/060012v3).

Essentially, there are two things you need to perform preranked gene set enrichment analysis with fgsea.

The fgsea package takes in two main arguments.

The first one is pathways: a list of gene sets or pathways to check

The second one is stats: a named vector of our genes of interest we want to perform GSEA on. The gene names must be the same as the ones in pathways! So make sure you are using the same nomenclature (i.e., gene IDs or Ensemble IDs).

That’s it!

There are additional parameters you may be interested in such as:

minSize: minimal size of a gene set to test (all pathways below the threshold are excluded)
maxSize: the same, but a maximum threshold.
scoreType: the default is ‘std’, where the enrichment score is computed normally, but you can also use one-tailed tests if you are only interested in positive (‘pos’) enrichment – so only pathways that are overrepresented – or negative (‘neg’) enrichment – to only get pathways that are underrepresented in your list of genes.

```R
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("fgsea")

library(fgsea)
library(ggplot2)
library(data.table)

# Example input data: ranked gene list
ranked_genes <- fread("ranked_genes.csv")$ranked_score

# Load gene sets (e.g., from MSigDB or KEGG)
pathways <- gmtPathways("pathways.gmt")

# Run fgsea
fgsea_results <- fgsea(pathways = pathways, stats = ranked_genes, nperm = 1000)

# Visualize results
ggplot(fgsea_results, aes(reorder(pathway, NES), NES)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "GSEA Results", x = "Pathway", y = "Normalized Enrichment Score")
```
> [!IMPORTANT]
> **How to Interpret GSEA Results**
> - *NES (Normalized Enrichment Score)*: Higher scores indicate more significant enrichment.
> - *p-value*: Determines the statistical significance of the enrichment.
> - *FDR (False Discovery Rate)*: A corrected p-value to account for multiple testing.
### 2.3 No-Code Alternative: GOrilla
GOrilla (Gene Ontology enRIchment anaLysis and visuaLizAtion tool) enables enrichment analysis through a web interface.
1. Access [GOrilla](http://cbl-gorilla.cs.technion.ac.il/).
2. Upload a ranked gene list.
3. Select the organism and analysis type (e.g., single-ranked list).
4. Explore enriched GO terms and export results.

---

## 3. Protein-Protein Interaction (PPI) Networks

Protein-Protein Interaction (PPI) networks represent the interactions between proteins, providing insights into biological processes and pathways. By integrating your list of differentially expressed genes (DEGs) into a PPI network, you can identify functional relationships and discover novel biological insights.

#### Steps for Building PPI Networks

1. **Identify the DEGs**: Start with a list of genes with significant differential expression.
2. **Obtain interaction data**: Use databases like STRING, BioGRID, or Pathway Commons to get protein interaction information.
3. **Construct the network**: Use R or external software to build and visualize the network.

### 3.1 Building PPI Networks in R with `STRINGdb`
```R
# Install and load STRINGdb for PPI analysis
BiocManager::install("STRINGdb")
library(STRINGdb)

# Initialize STRINGdb
string_db <- STRINGdb$new(species = 9606)  # Homo sapiens

# Map gene IDs and retrieve interactions
genes <- c("BRCA1", "TP53", "EGFR")  # Example genes
mapped_genes <- string_db$map(genes, "gene")
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

# Visualize in R or export to Cytoscape
```

### 3.2 Visualizing PPI Networks with Cytoscape
1. Download [Cytoscape](https://cytoscape.org/).
2. Import STRINGdb interaction file.
3. Use built-in apps (e.g., stringApp) to enhance visualizations.

### 3.3 No-Code Alternative Tools
- [STRING Database](https://string-db.org/): Explore interactions directly on the web interface.

---

## 4. Drug Repurposing
### 4.1 Introduction
Drug repurposing identifies existing drugs that could target pathways or genes of interest. RNA-seq data helps uncover potential targets by linking DEGs to known drug-gene interactions.

### 4.2 Drug Repurposing in R
```R
# Install and use drugfindR
devtools::install_github("MaayanLab/drugfindR")
library(drugfindR)

# Example: Search for drugs targeting upregulated genes
upregulated_genes <- c("EGFR", "VEGFA", "STAT3")
drug_results <- drugfindR::find_drugs(upregulated_genes)

# View results
head(drug_results)
```

### 4.3 No-Code Alternative: LINCS L1000 Platform
1. Access [LINCS L1000](https://clue.io/).
2. Upload a list of DEGs or ranked genes.
3. Discover drug perturbations that match or reverse your signature.

---

## 5. Biological Interpretation and Reporting
### 5.1 Strategies for Interpretation
- Focus on pathways and networks relevant to your research question.
- Validate findings with additional datasets or literature.

### 5.2 Communicating Results
- Use clear visualizations (e.g., enrichment plots, PPI diagrams).
- Report findings with biological context, emphasizing their relevance to disease mechanisms or therapeutic opportunities.

---

## Additional Resources
- [Pathway Enrichment Analysis (PEA)](https://biostatsquid.com/pathway-enrichment-analysis-explained/)
-  [Gene Set Enrichment Analysis](https://biostatsquid.com/gene-set-enrichment-analysis/)
- [fgsea Documentation](https://bioconductor.org/packages/release/bioc/html/fgsea.html)
- [Using fgsea package](https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html)
- [Proteomics Data Analysis in R](https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/pathway-analysis.html)
- [GOrilla Tool](http://cbl-gorilla.cs.technion.ac.il/)
- [STRING Database](https://string-db.org/)
- [LINCS L1000 Platform](https://clue.io/)

---

## References

- Subramanian, A., et al. (2005). Gene Set Enrichment Analysis: A Knowledge-Based Approach for Interpreting Genome-Wide Expression Profiles. *Proceedings of the National Academy of Sciences*.
- Szklarczyk, D., et al. (2019). STRING v11: protein–protein association networks with increased coverage, supporting functional discovery in genome-wide experimental datasets. *Nucleic Acids Research*.
- Gene Ontology Consortium. (2019). The Gene Ontology Resource: 20 years and still GOing strong. *Nucleic Acids Research*.
- Eran Eden*, Roy Navon*, Israel Steinfeld, Doron Lipson and Zohar Yakhini. "GOrilla: A Tool For Discovery And Visualization of Enriched GO Terms in Ranked Gene Lists", BMC Bioinformatics 2009, 10:48.
- Eran Eden, Doron Lipson, Sivan Yogev, Zohar Yakhini. "Discovering Motifs in Ranked Lists of DNA sequences", PLoS Computational Biology, 3(3):e39, 2007.
