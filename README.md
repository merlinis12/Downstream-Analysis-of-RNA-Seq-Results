
# Downstream Analysis of RNA-Seq Results in R: GSEA, PPI Networks, and Biological Interpretation

## Workshop Overview
After identifying differentially expressed genes as we did in the [previous workshop]((https://github.com/merlinis12/RNA-Seq-Data-Analysis-in-R/wiki)), it’s critical to interpret these findings in a biological context. 
This workshop explores downstream analyses you can perform with differential gene expression analysis results in R.
You’ll learn how to:
- Conduct [Gene Set Enrichment Analysis (GSEA)](#2-gene-set-enrichment-analysis-gsea) for pathway-level insights.
- Build [Protein-Protein Interaction (PPI) networks](#3-protein-protein-interaction-ppi-networks) to identify functional connections between genes.
- Explore drug repurposing opportunities based on RNA-seq data.
Hands-on exercises will focus on tools like `fgsea` in R and software such as GOrilla and Cytoscape.

---

## 1. Introduction to Downstream Analysis
Downstream analysis is essential for understanding the biological significance of differentially expressed genes (DEGs). 
We aim to connect RNA-seq data to larger biological systems, providing insights into gene function, pathways, and therapeutic opportunities.

---

## 2. Gene Set Enrichment Analysis (GSEA)

Comment on the ranked list of degs.

Comment on the gene id conversion.

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
