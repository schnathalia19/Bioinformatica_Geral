# Tarefa 2 — Análise de Bulk RNA-seq

## 1. Dataset

| Recurso | Identificador |
|---|---|
| GEO | [GSE348127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348127) |
| SRA Run Selector | [PRJNA1532893](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1532893&o=acc_s%3Aa) |

## 2. Alinhamento e contagem

### 2.1. Download dos arquivos FASTQ

```bash
mkdir -p fastq

for r in SRR40776939 SRR40776940 SRR40776941 SRR40776942 SRR40776943 SRR40776944; do
  echo "Baixando $r..."

  fastq-dump -X 5000000 \
    --split-3 \
    --gzip \
    --skip-technical \
    --readids \
    -O fastq "$r"
done
```

> A opção `-X 5000000` limita o download a 5 milhões de bases por amostra.

### 2.2. Controle de qualidade

```bash
mkdir -p qc

fastqc fastq/*.fastq.gz -o qc
multiqc qc -o qc
```

### 2.3. Genoma de referência e anotação

```bash
wget [https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz](https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz)
tar xzf grch38_genome.tar.gz

wget [https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz](https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz)
gunzip Homo_sapiens.GRCh38.112.gtf.gz
```

### 2.4. Alinhamento e processamento dos arquivos BAM

```bash
mkdir -p bam

for r in SRR40776939 SRR40776940 SRR40776941 SRR40776942 SRR40776943 SRR40776944; do
  hisat2 -p 8 \
    --rna-strandness RF \
    -x grch38/genome \
    -1 fastq/${r}_1.fastq.gz \
    -2 fastq/${r}_2.fastq.gz \
    --summary-file bam/${r}.hisat2.txt \
  | samtools sort -@ 4 -o bam/${r}.bam

  samtools index bam/${r}.bam
  samtools flagstat bam/${r}.bam > bam/${r}.flagstat.txt
done

grep "overall alignment rate" bam/*.hisat2.txt
```

### 2.5. Contagem de leituras por gene

```bash
mkdir -p counts

featureCounts -p \
  --countReadPairs \
  -s 2 \
  -a Homo_sapiens.GRCh38.112.gtf \
  -o counts/counts.txt \
  bam/*.bam

cat counts/counts.txt.summary
```

## 3. Expressão diferencial e GSEA em R

> **Antes de executar:** crie o objeto `dds` a partir da matriz de contagens e dos metadados das amostras. Os metadados precisam identificar a condição (`condition`) e o doador (`donor`) de cada amostra. Essa etapa não está incluída no código abaixo.

### 3.1. Carregar os pacotes

```r
library(tidyverse)
library(DESeq2)
library(apeglm)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
```

### 3.2. Análise com DESeq2

```r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "BaP", "CTRL")
)

res_shrunk <- lfcShrink(
  dds,
  coef = "condition_BaP_vs_CTRL",
  type = "apeglm"
)

res_shrunk$symbol <- mapIds(
  org.Hs.eg.db,
  keys = row.names(res_shrunk),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
```

### 3.3. PCA

```r
vsd <- vst(dds, blind = FALSE)

pca_plot <- plotPCA(
  vsd,
  intgroup = c("condition", "donor")
) +
  ggtitle("PCA - Efeito do Tratamento BaP") +
  theme_minimal()

print(pca_plot)
```

### 3.4. Volcano plot

```r
volcano_plot <- EnhancedVolcano(
  res_shrunk,
  lab = res_shrunk$symbol,
  x = "log2FoldChange",
  y = "padj",
  title = "Volcano Plot: BaP vs CTRL",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0
)

print(volcano_plot)
```

### 3.5. Heatmap dos 40 principais genes

```r
top_genes <- head(
  order(res_shrunk$padj, decreasing = FALSE),
  40
)

mat <- assay(vsd)[top_genes, ]
mat <- mat - rowMeans(mat)
rownames(mat) <- res_shrunk$symbol[top_genes]

pheatmap(
  mat,
  annotation_col = as.data.frame(
    colData(dds)[, "condition", drop = FALSE]
  ),
  main = "Top 40 DEGs (Heatmap)",
  fontsize_row = 8
)
```

### 3.6. GSEA com vias Hallmark

```r
res_gsea <- as.data.frame(res_shrunk) %>%
  filter(!is.na(symbol) & !is.na(log2FoldChange))

gene_list <- res_gsea$log2FoldChange
names(gene_list) <- res_gsea$symbol
gene_list <- sort(gene_list, decreasing = TRUE)

m_t2g <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(
  gene_list,
  TERM2GENE = m_t2g,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)
```

### 3.7. Dotplot das vias enriquecidas

```r
dot_plot <- dotplot(
  gsea_res,
  split = ".sign",
  showCategory = 10
) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA Dotplot - Vias Hallmark")

print(dot_plot)
```

### 3.8. Gráfico das três principais vias

```r
gsea_top3 <- gseaplot2(
  gsea_res,
  geneSetID = 1:3,
  title = "Top 3 Vias Enriquecidas"
)

print(gsea_top3)
```
