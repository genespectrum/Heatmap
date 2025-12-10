### Load libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(biomaRt)

### 1️⃣ Load count file
file_path <- file.choose()
counts <- read.delim(file_path, header = TRUE, row.names = 1)

### 2️⃣ Remove the 'Length' column
counts <- counts[, !colnames(counts) %in% "Length"]

### 3️⃣ Select samples for comparison
samples <- c("Control.1","Control.2","Control.3",
             "Cys247X.1","Cys247X.2","Cys247X.3")
counts <- counts[, samples]

### 4️⃣ Create sample metadata
condition <- factor(c(rep("Control",3), rep("Cys247X",3)))
coldata <- data.frame(row.names = samples, condition)

### 5️⃣ Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)

### 6️⃣ Extract DE results
res <- results(dds)

### 7️⃣ Filter significant DE genes
sig_res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
sig_res <- sig_res[order(sig_res$padj), ]

### 8️⃣ Select top significant genes (limit to top 100)
top_genes <- rownames(sig_res)[1:min(30, nrow(sig_res))]

### 9️⃣ Variance Stabilizing Transformation
vsd <- vst(dds, blind = FALSE)

### 🔟 Gene ID → Gene Name mapping using biomaRt
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

gene_ids <- rownames(vsd)

gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = gene_ids,
                  mart = mart)

### 1️⃣1️⃣ Ensure uniqueness and match order
gene_map <- gene_map %>% distinct(ensembl_gene_id, .keep_all = TRUE)

mapped_symbols <- gene_map$hgnc_symbol[
  match(gene_ids, gene_map$ensembl_gene_id)
]

### Replace blanks with Ensembl IDs
#mapped_symbols[mapped_symbols == ""] <- gene_ids[mapped_symbols == ""]
# Replace NA or blank gene names with Ensembl IDs
idx <- is.na(mapped_symbols) | mapped_symbols == ""
mapped_symbols[idx] <- gene_ids[idx]

### Assign new rownames to the VST matrix
rownames(vsd) <- mapped_symbols

### Convert top_genes (original IDs) → gene names
top_gene_names <- mapped_symbols[match(top_genes, gene_ids)]

### 1️⃣2️⃣ Extract VST matrix for top significant genes
mat <- assay(vsd)[top_gene_names, ]

### 1️⃣3️⃣ Plot heatmap with gene names
pheatmap(mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 10,
         main = "Heatmap of Significant DE Genes (padj < 0.05, Gene Names)")

