# <span style="color: LightSkyBlue">Downstream Analysis of RNA-Seq Results in R: GSEA, PPI Networks, and Biological Interpretation</span>

## <span style="color: greenyellow">Workshop Overview</span>

After identifying differentially expressed genes (DEGs) as we did in the [previous workshop](https://github.com/merlinis12/RNA-Seq-Data-Analysis-in-R/wiki), it’s critical to interpret these findings in a biological context. In this workshop, we will explore downstream analyses you can perform using RNA-seq data.

### Key Learning Goals:
- Conduct [Gene Set Enrichment Analysis (GSEA)](#2-gene-set-enrichment-analysis-gsea) for pathway-level insights.
- Build [Protein-Protein Interaction (PPI) Networks](#3-protein-protein-interaction-ppi-networks) to explore functional connections between genes.
- Explore drug repurposing opportunities based on RNA-seq data.

---

## <span style="color: LightSkyBlue">1. Introduction to Downstream Analysis</span>

Differential gene expression (DGE) analysis is an essential step in RNA-seq data analysis. The goal is to identify DEGs between two conditions. For example, in the [previous workshop](https://github.com/merlinis12/RNA-Seq-Data-Analysis-in-R), we studied the difference in gene expression between airway smooth muscle cells from healthy individuals and those treated with dexamethasone.

### The Challenge:
After performing DGE, you might be left with a long list of DEGs.

:pinched_fingers: How do we start interpreting these results?

:pinched_fingers: Is there a way to simplify and interpret hundreds of DEGs at once?

> :mechanical_arm: **Solution:**  
>
> One common approach is **Pathway Enrichment Analysis (PEA)**, which summarizes large gene lists into a more easily interpretable list of biological pathways. Instead of reviewing thousands of genes, you’ll end up with a shorter list of relevant pathways.

---

## <span style="color: LightSkyBlue">2. Pathway Enrichment Analysis (PEA)</span>

> :brain: TO REMEMBER
>
> **Pathway enrichment analysis** compares your gene list to a background list to check if certain pathways are overrepresented.

### PEA Workflow:
1. **Gene List**: This could be your differentially expressed genes which you want to summarize.
2. **Background Gene List**: For example, all genes in the human genome.
3. **Gene Sets**: Groups of related genes associated with specific pathways or diseases.

You will learn more about each component as we go along.

**INSERT THE PICTURE!!!!!!!!!!!!!!!!!!!!!!!!**

PEA helps answer questions like:
> **Is our list of differentially expressed genes enriched with genes involved in the IL-6 synthesis pathway?**

To test this, we use a contingency table comparing the fraction of genes in the pathway vs. those outside of it (background).

Here’s an example using 30 genes:
- 15 DEGs
- 12 out of those 15 DEGs are involved in IL-6 production. We could quite confidently say that our gene list is enriched with genes involved in IL-6 production.

We can statistically test this using Fisher’s Exact Test to check if the IL-6 production pathway is overrepresented.

>  [!IMPORTANT]  
> A low p-value indicates enrichment, but **check the specific genes** involved—are they positive or negative regulators?

#### Pathway Enrichment in Action:
> [!CAUTION]  
> A low p-value **doesn’t mean** the pathway is upregulated. You need to check if the genes are promoting or inhibiting the pathway.

---

### Multiple Testing Correction

PEA involves testing multiple pathways, which introduces a risk of finding significant results by chance. This is addressed by using a **multiple testing correction** method, such as **Benjamini-Hochberg (BH)**.

> [!NOTE]  
> **Summary**:  
> PEA compares your long gene list to a background set and determines if certain pathways are overrepresented. It involves checking the fraction of genes annotated to a specific gene set (like Gene Ontology (GO) terms) and comparing it to the genome-wide distribution.



---

## <span style="color: LightSkyBlue">3. Gene Sets and Databases</span>

### Common Gene Sets Databases:
- **Gene Ontology** [(GO)](https://geneontology.org/): Focuses on biological processes, molecular functions, and cellular components (many species including humans).
- Kyoto Encyclopedia of Genes and Genomes ([KEGG](https://www.genome.jp/kegg/)): Primarily concerned with metabolic pathways (many species including humans).
- [**Reactome**](https://reactome.org/): A curated database of human molecular pathways (only humans).
  

| **Aspect**                | **Gene Ontology (GO)**                               | **KEGG**                                         | **Reactome**                                        |
|---------------------------|-----------------------------------------------------|-------------------------------------------------|-----------------------------------------------------|
| **Focus Area**             | Biological processes, molecular functions, cellular components | Metabolic pathways and molecular networks       | Curated human molecular pathways                   |
| **Scope**                  | Universal; applies to all organisms                 | Primarily focused on metabolic pathways         | Focuses on human biological pathways                |
| **Categories**             | Biological Processes (BP), Molecular Functions (MF), Cellular Components (CC) | Metabolic, signaling pathways, disease pathways | Molecular interactions, cellular processes, pathways|
| **Usage**                  | Broadly used for gene functional annotations       | Primarily used for pathway mapping in metabolism | Used for pathway visualization and analysis         |
| **Data Type**              | Terms representing biological concepts              | Pathways, networks, and diseases                | Pathways and interactions                           |
| **Example Terms**          | "Apoptosis", "Cell cycle", "DNA binding"            | "Glycolysis", "Cytokine signaling", "Cancer"    | "Signal transduction", "Gene expression regulation" |
| **Application**            | Gene annotation and enrichment analysis             | Metabolic profiling and drug development        | Disease research, drug discovery, and systems biology|
| **Method of Representation** | Hierarchical terms with relationships             | Pathways mapped with enzymes, metabolites, and genes | Curated pathway maps with gene interactions         |
| **Available Tools**        | GO-Term Finder, GOseq, clusterProfiler              | KEGG Mapper, Pathview                           | Reactome Pathway Browser, Pathway Analysis          |
| **Human-Specific**         | No, applicable to all organisms                     | No, but it includes many human pathways         | Yes, Reactome is focused solely on human data       |
| **URL**                    | [Gene Ontology](https://www.geneontology.org/)       | [KEGG](https://www.kegg.jp/)                     | [Reactome](https://reactome.org/)                   |


#### Gene Ontology (GO)
Gene Ontology is a controlled vocabulary for describing gene functions. It consists of three ontologies:
1. **Biological Processes (BP)**: Describes large biological pathways or events (e.g., "cell cycle").
2. **Molecular Functions (MF)**: Describes molecular activities (e.g., "protein kinase activity").
3. **Cellular Components (CC)**: Describes the location within the cell (e.g., "plasma membrane").

> [!NOTE]  
> **GO Terms Hierarchy**: GO terms are organized hierarchically, where broader terms are parents to more specific ones.

### Selecting Background Genes
Choosing the right background gene set is crucial for accurate pathway enrichment analysis. There are gene sets for diseases, tissues, transcription factors, and more.

> [!WARNING]  
> Carefully select your background genes based on your experimental setup to ensure accurate results.

---



## <span style="color: LightSkyBlue"> 4. Gene List</span>

Gene expression analysis results often look like this:

- Some genes are **downregulated**.
- Some are **upregulated**.
- Some **don’t change** at all.

Not all changes are significant, and some might not be differentially expressed at all. 

There are some methods, called **Over-Representation Analysis (ORA) methods** that allow to first filter your results, by significance and fold change, to keep only differentially expressed genes.


The significance of each pathway is measured by calculating the probability that the observed number of DE genes in a given pathway occurred by chance. Lower p-values indicate overrepresentation.

> [!CAUTION]  
> **Overrepresentation analysis methods** heavily depend on the cutoffs used to select DE genes.



:pinched_fingers: **Is there a more objective or unbiased way to handle this?**

> :mechanical_arm: **Solution:**  
>
> **Gene Set Enrichment Analysis (GSEA)**, a great method for pathway enrichment that takes both significance and fold change into account.





---

## <span style="color: aqua;">5. Gene Set Enrichment Analysis (GSEA)</span>

Gene Set Enrichment Analysis (GSEA) is a method used to determine whether a predefined set of genes shows statistically significant differences in expression between two biological conditions. GSEA does not rely on a threshold (like fold-change) and instead assesses whether the genes in a set are overrepresented at the extremes of a ranked list of genes.

**Steps in GSEA**

- **Rank the genes** based on their differential expression.
- **Predefine gene sets**: These could be biological pathways, gene ontology (GO) terms, or other gene collections.
- **Calculate enrichment scores** for each gene set.
- **Assess statistical significance** using permutation testing.


> [!IMPORTANT]  
> The steps to perform **Gene Set Enrichment Analysis** are very similar to those of overrepresentation analysis methods. The key difference is that GSEA takes a **ranked list of genes** as input rather than a filtered list.

:pinched_fingers: **How does it work?**

Genes are ranked by some score that combines **differential expression** and **direction of change**. A common ranking method is:

$$**Ranking = Sign(log2FoldChange) * -log10(p-value)**$$

This gives us a ranked list that orders genes by significance and the direction of change (upregulated or downregulated).

The results of GSEA are similar to overrepresentation analysis, but with the added benefit that **genes are not filtered** beforehand, and it considers the magnitude and significance of the changes.

> [!NOTE]  
> You can just use the $log2(fold-change)$ as a ranking – but it will not take into account genes with a large fold change but non-significant. However, if you have already selected your significant genes, this may be a good option for you.


**Key Benefits of GSEA:**

- GSEA allows you to **rank genes by significance and fold change**.
- It does not require filtering genes beforehand.
- **GSEA compares the ranked list** with a list of background genes and determines which pathways are overrepresented based on the ranked list.


> [!NOTE]  
> There is a newer generation of pathway enrichment methods called **Topology-based (TB) methods**, which take into account gene dependencies and interactions (more [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3146-1)).


### 5.1 GSEA in R with `fgsea`

The `fgsea` package in R provides a fast and efficient method for performing preranked Gene Set Enrichment Analysis (GSEA). It enables precise and rapid calculation of extremely low GSEA p-values across a collection of gene sets. The p-value estimation employs an adaptive multi-level split Monte Carlo approach. For a detailed explanation of the algorithm, refer to [this preprint](https://www.biorxiv.org/content/10.1101/060012v3).

Essentially, there are two things you need to perform preranked gene set enrichment analysis with `fgsea`.

The `fgsea` package takes in two main arguments:

1. **pathways**: A list of gene sets or pathways to check.
2. **stats**: A named vector of your genes of interest for GSEA. The gene names must match the ones in the pathways list. Make sure you are using the same nomenclature (e.g., gene IDs or Ensembl IDs).

That’s it!

There are additional parameters you may be interested in, such as:

- **minSize**: The minimal size of a gene set to test. Pathways below this threshold are excluded.
- **maxSize**: The maximum size of a gene set to test. Pathways above this threshold are excluded.
- **scoreType**: The default is ‘std’, where the enrichment score is computed normally. However, you can also use one-tailed tests if you are only interested in positive (‘pos’) enrichment – so only pathways that are overrepresented – or negative (‘neg’) enrichment – to only get pathways that are underrepresented in your list of genes.



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
> - 
### 5.2 No-Code Alternative: GOrilla
GOrilla (Gene Ontology enRIchment anaLysis and visuaLizAtion tool) enables enrichment analysis through a web interface.
1. Access [GOrilla](http://cbl-gorilla.cs.technion.ac.il/).
2. Upload a ranked gene list.
3. Select the organism and analysis type (e.g., single-ranked list).
4. Explore enriched GO terms and export results.

---

## <span style="color: aqua;">6. Protein-Protein Interaction (PPI) Networks</span>

## <span style="color: aqua;">6. Protein-Protein Interaction (PPI) Networks</span>

Protein-Protein Interaction (PPI) networks are powerful tools for understanding how proteins work together within cells. By integrating **differentially expressed genes (DEGs)** into these networks, we can uncover functional relationships and gain deeper insights into biological processes and pathways.

:pinched_fingers: **What is STRING?**

**STRING** (Search Tool for the Retrieval of Interacting Genes/Proteins) is a comprehensive database that provides information about both **functional** and **physical protein-protein interactions** across different species. It combines various types of evidence, including:

- **Experimental data**: Direct interactions observed in laboratory experiments.
- **Computational predictions**: Based on genomic context and co-expression.
- **Text mining**: Analyzing scientific literature to predict interactions.
- **Pathway annotations**: Interactions tied to specific biological pathways.

STRING also allows you to perform **clustering** and **pathway enrichment analysis**, and it visualizes your protein interaction data to help identify key biological processes.


Ready to get started? Let’s dive into some hands-on PPI analysis using the [STRING Database](https://string-db.org/).



---

## <span style="color: aqua">7. Drug Repurposing<span>
### 4.1 Introduction
Drug repurposing identifies existing drugs that could target pathways or genes of interest. RNA-seq data helps uncover potential targets by linking DEGs to known drug-gene interactions.

To perform drug repurposing in R using a list of differentially expressed genes, you can leverage packages like "DrugVsDisease" or "enrichR" to compare your gene expression signature against large drug perturbation datasets like the [LINCS L1000](https://lincsportal.ccs.miami.edu/signatures/datasets/LDG-1188), identifying drugs that show a significantly reversed expression pattern compared to your disease gene set, essentially indicating potential therapeutic applications for those drugs in the disease context. 

**Key steps:**

**Prepare your gene list:**
- Ensure your differentially expressed genes are ranked based on their fold change or significance level (p-value).
- Consider filtering for high-confidence genes based on adjusted p-values.

**Access a drug perturbation database:**
- Utilize public datasets like the LINCS L1000, which provides gene expression profiles for a wide range of drugs across different cell lines. 

**Calculate similarity scores:**
- Use a method like the Kolmogorov-Smirnov test to compare the ranked gene expression profile of your disease signature to each drug's gene expression profile from the database. 
- A high similarity score indicates that the drug might reverse the disease-related gene expression pattern. 

**Enrichment analysis:**
- Employ packages like "enrichR" to identify pathways or gene sets significantly enriched within your differentially expressed genes, which can further guide drug repurposing by identifying potential therapeutic targets. 


**R packages to consider:**

- `DrugVsDisease`: A dedicated package for drug repurposing analysis using gene expression signatures, calculating enrichment scores based on drug-induced gene expression changes. 
- `clusterProfiler`: Offers gene set enrichment analysis functionalities that can be used to identify pathways associated with your differentially expressed genes. 

```R
# Load libraries

library(DrugVsDisease)

library(enrichR)



# Load your differentially expressed genes (DEG) data

DEG <- read.csv("DEG_list.csv") 



# Access LINCS L1000 data (you might need to download and process this separately)

lincs_data <- load_lincs_data() 



# Calculate similarity scores for each drug against your DEG list

drug_scores <- calculate_similarity(DEG, lincs_data)



# Identify top potential repurposing drugs

top_drugs <- drug_scores %>% 

    arrange(similarity_score) %>% 

    head(n = 10) 



# Perform pathway enrichment analysis on your DEG list 

enrichment_results <- enrichPathway(DEG$gene_id,  

                                   organism = "Homo sapiens", 

                                   pvalueCutoff = 0.05) 



# Analyze the results to identify potential therapeutic pathways 

plot_enrichment(enrichment_results)

```

> [!IMPORTANT]
> - Interpret results based on the known biology of your disease and the target drugs.
> - Further experimental validation is necessary to confirm potential drug repurposing candidates identified computationally.

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
## <span style="color: aqua">6. VALIDATION, VALIDATION, VALIDATION</span>

---

## Widely Used Bioinformatics Tools for Biological Pathway Analysis
- **Enrichr**: A tool for gene set enrichment analysis and visualization. [Enrichr](https://maayanlab.cloud/Enrichr/)
- **DAVID (Database for Annotation, Visualization, and Integrated Discovery)**: A comprehensive set of functional annotation tools. [DAVID](https://david.ncifcrf.gov/)
- **Metascape**: A versatile platform for gene annotation and analysis. [Metascape](https://metascape.org/gp/index.html)
- **GSEA (Gene Set Enrichment Analysis)**: A method for identifying enriched gene sets. [GSEA](http://www.gsea-msigdb.org/gsea/index.jsp)
- **GSVA (Gene Set Variation Analysis)**: A method for analyzing gene set variation. [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)
- **Cytoscape**: A tool for visualizing and analyzing biological networks. [Cytoscape](https://cytoscape.org/)


##  Resources
- [Pathway Enrichment Analysis (PEA)](https://biostatsquid.com/pathway-enrichment-analysis-explained/)
- [Gene Set Enrichment Analysis](https://biostatsquid.com/gene-set-enrichment-analysis/)
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
- Rosati D, Palmieri M, Brunelli G, Morrione A, Iannelli F, Frullanti E, Giordano A. Differential gene expression analysis pipelines and bioinformatic tools for the identification of specific biomarkers: A review. Comput Struct Biotechnol J. 2024 Mar 1;23:1154-1168. doi: 10.1016/j.csbj.2024.02.018. PMID: 38510977; PMCID: PMC10951429.
