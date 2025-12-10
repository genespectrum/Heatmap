# Heatmap of Significant Differentially Expressed Genes
### This repository provides a simple and clear guide to help you generate a publication-ready heatmap of significant differentially expressed (DE) genes using DESeq2, biomaRt, and pheatmap in R.

## Overview
This workflow identifies significant DE genes between two conditions, transforms the data for better visualization, converts Ensembl gene IDs to gene names, and generates a heatmap showing expression patterns. This type of heatmap is commonly used in transcriptomics studies to highlight clusters of genes that behave similarly across samples.
The final heatmap displays:
- Variance-stabilized expression values
- Gene names instead of Ensembl IDs
- Clustering across genes and samples
- Clear separation between experimental conditions
## Dataset Requirements
You need a raw read count file containing Ensembl gene IDs as row names and sample names as column names. The file must not contain any non-numeric columns except optional metadata like gene length. Make sure the sample names match exactly with the names you want to analyze.
Your dataset must include:
- A column for each sample
- Rows representing genes
- Count values (integers)
## Step-by-Step Tutorial

1. Install the Required Libraries
   
   You begin by Installing the necessary Packages such as DESeq2, pheatmap, dplyr, ggplot2 and biomaRt. These packages help you perform differential expression analysis, transform gene expression values, annotate gene IDs, and generate the heatmap.

2. Upload the Raw Count File
   
   Select your count file using a file chooser window, and R reads it into a data frame. The file must have gene IDs in the first column and sample names in the remaining columns.

3. Remove the Length Column

   If your file contains a column named “Length,” it is removed to ensure that only raw counts remain. This step prevents errors during DESeq2 processing.

4. Select Samples for Comparison

   You specify which samples belong to your control and treatment groups. Only these samples are kept for downstream analysis.

5. Create the Sample Metadata

   A simple metadata table is created to indicate the condition of each sample. This table is used by DESeq2 to understand the experimental design.

6. Run the DESeq2 Pipeline

   DESeq2 processes the count data, normalizes it, estimates dispersion, and identifies genes that are significantly different between the two groups.

7. Extract Significant Results

   The DE result table is filtered to keep only genes with strong statistical significance and meaningful fold changes. The significant genes are then sorted by adjusted p-value.

8. Select the Top Significant Genes

   You choose the top genes based on significance. These genes will be used in the heatmap for clearer visualization.

9. Apply Variance Stabilizing Transformation (VST)

   VST is applied to reduce noise and make sample-to-sample differences more comparable. This helps the heatmap show clearer clustering.

10. Map Gene IDs to Gene Names

    Using biomaRt, each Ensembl gene ID is mapped to its corresponding gene symbol. If a symbol is missing, the original ID is kept.

11. Match Gene Names Correctly

    The final list of gene names is matched with the transformed matrix so that the row names in the heatmap are properly labeled.

12. Extract the Expression Matrix

    The variance-stabilized expression values for the top genes are extracted. This matrix forms the input for the heatmap.

13. Generate the Heatmap

    A heatmap is created using pheatmap. Rows and columns are clustered, values are scaled, and a clear title is added to highlight significant DE genes.

## Use Cases

1. RNA-seq Differential Expression Analysis
   Heatmaps help visualize expression patterns of top significant genes between two conditions.

2. Biomarker Discovery

   Researchers use heatmaps to identify groups of genes that may act as biomarkers based on clustering patterns.

3. Pathway-Level Expression Patterns

   Heatmaps show how genes associated with the same pathway behave across samples.

4. Sample Clustering and Quality Assessment

   The heatmap reveals how similar or different samples are, helping identify batch effects or outliers.
