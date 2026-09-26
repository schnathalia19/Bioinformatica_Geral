# Tarefa 2 — Bulk RNA-seq

## 1. Escolha do dataset
- GEO: [GSE348127](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348127)
- SRA Run Selector: [PRJNA1532893](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1532893&o=acc_s%3Aa)

## 2. Alinhamento e contagem

### Download dos FASTQ
```bash
mkdir -p fastq
for r in SRR40776939 SRR40776940 SRR40776941 SRR40776942 SRR40776943 SRR40776944; do
  echo "Baixando $r..."
  fastq-dump -X 5000000 --split-3 --gzip --skip-technical --readids -O fastq "$r"
done
```

### Controle de qualidade
```bash
mkdir -p qc
fastqc fastq/*.fastq.gz -o qc
multiqc qc -o qc
```

### Índice do genoma e anotação
```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xzf grch38_genome.tar.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz
```

### Alinhamento + sort + index + flagstat
```bash
mkdir -p bam
for r in SRR40776939 SRR40776940 SRR40776941 SRR40776942 SRR40776943 SRR40776944; do
  hisat2 -p 8 --rna-strandness RF -x grch38/genome \
    -1 fastq/${r}_1.fastq.gz -2 fastq/${r}_2.fastq.gz \
    --summary-file bam/${r}.hisat2.txt \
  | samtools sort -@ 4 -o bam/${r}.bam
  samtools index bam/${r}.bam
  samtools flagstat bam/${r}.bam > bam/${r}.flagstat.txt
done
grep "overall alignment rate" bam/*.hisat2.txt
```

### Contagem
```bash
mkdir -p counts
featureCounts -p --countReadPairs -s 2 \
  -a Homo_sapiens.GRCh38.112.gtf \
  -o counts/counts.txt bam/*.bam
cat counts/counts.txt.summary
```

## 3. Expressão diferencial e GSEA em R

### Pacotes
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "EnhancedVolcano", "org.Hs.eg.db",
                       "clusterProfiler", "enrichplot", "edgeR", "limma"))
install.packages(c("tidyverse", "pheatmap", "msigdbr"))
```

### Importar contagens e montar o colData
```r
library(tidyverse); library(DESeq2); library(apeglm); library(pheatmap)
library(EnhancedVolcano); library(org.Hs.eg.db); library(clusterProfiler)
library(enrichplot); library(msigdbr)

fc <- read.delim("counts/counts.txt", comment.char = "#", check.names = FALSE)
counts <- as.matrix(fc[, 7:ncol(fc)])
rownames(counts) <- fc$Geneid
colnames(counts) <- sub("\\.bam$", "", basename(colnames(counts)))

coldata <- data.frame(
  run       = c("SRR40776944", "SRR40776943", "SRR40776942",
                "SRR40776941", "SRR40776940", "SRR40776939"),
  condition = factor(rep(c("CTRL", "BaP"), each = 3), levels = c("CTRL", "BaP")),
  donor     = factor(rep(c("NK187", "NK110", "12SKera006"), 2))
)
rownames(coldata) <- coldata$run

# IMPORTANTE: mesma ordem na matriz e no colData
counts <- counts[, rownames(coldata)]
stopifnot(all(colnames(counts) == rownames(coldata)))

# Matriz de contagens brutas (entregável)
write.csv(counts, "matriz_contagens_brutas.csv")
```

### DESeq2
```r
# Desenho pareado: controla o efeito da doadora
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ donor + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]   # filtro de baixa expressão
dds <- DESeq(dds)

res    <- results(dds, contrast = c("condition", "BaP", "CTRL"), alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "condition_BaP_vs_CTRL", type = "apeglm")
summary(res)

res_df <- as.data.frame(resLFC) |>
  rownames_to_column("ensembl") |>
  mutate(symbol = mapIds(org.Hs.eg.db, ensembl, "SYMBOL", "ENSEMBL", multiVals = "first"),
         stat   = res[ensembl, "stat"]) |>
  arrange(padj)

degs <- filter(res_df, padj < 0.05, abs(log2FoldChange) >= 1)
nrow(degs)

write.csv(res_df, "resultados_DESeq2.csv", row.names = FALSE)
```

### PCA
```r
vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("condition", "donor")) +
  geom_point(aes(shape = donor), size = 4) + theme_bw()
ggsave("PCA.png", width = 6, height = 4.5)

# PCA removendo o efeito da doadora (só para visualização)
vsd_bc <- vsd
assay(vsd_bc) <- limma::removeBatchEffect(assay(vsd), batch = vsd$donor,
                                          design = model.matrix(~ condition, colData(vsd)))
plotPCA(vsd_bc, intgroup = "condition") + theme_bw()
ggsave("PCA_sem_efeito_doadora.png", width = 6, height = 4.5)
```

### Volcano plot
```r
EnhancedVolcano(res_df, lab = res_df$symbol,
                x = "log2FoldChange", y = "padj",
                pCutoff = 0.05, FCcutoff = 1,
                title = "BaP vs DMSO", subtitle = "DESeq2 (apeglm)",
                ylab = bquote(~-Log[10] ~ italic(P)[adj]))
ggsave("volcano.png", width = 8, height = 7)
```

### Heatmap dos principais DEGs
```r
top <- head(filter(res_df, padj < 0.05), 50)
mat <- assay(vsd)[top$ensembl, ]
rownames(mat) <- ifelse(is.na(top$symbol), top$ensembl, top$symbol)

pheatmap(mat, scale = "row",
         annotation_col = coldata[, c("condition", "donor")],
         show_colnames = FALSE, fontsize_row = 7,
         filename = "heatmap_top50.png", width = 6, height = 9)
```

### GSEA (todos os genes, ranqueados pela estatística de Wald)
```r
ranks <- res_df |> filter(!is.na(stat)) |> select(ensembl, stat) |> deframe()
ranks <- sort(ranks, decreasing = TRUE)
set.seed(123)

# GO Biological Process (IDs Ensembl direto, sem conversão)
gse_go <- gseGO(ranks, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP",
                minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05,
                eps = 0, seed = TRUE)
gse_go <- setReadable(gse_go, org.Hs.eg.db, keyType = "ENSEMBL")

# Hallmarks do MSigDB
# (em versões antigas do msigdbr, use category = "H" em vez de collection = "H")
h <- msigdbr(species = "Homo sapiens", collection = "H") |>
  select(gs_name, ensembl_gene)
gse_h <- GSEA(ranks, TERM2GENE = h, pvalueCutoff = 0.05, eps = 0, seed = TRUE)

# Figuras
dotplot(gse_go, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
ggsave("GSEA_GO_dotplot.png", width = 10, height = 7)

dotplot(gse_h, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign)
ggsave("GSEA_Hallmark_dotplot.png", width = 10, height = 6)

gseaplot2(gse_h, geneSetID = 1:3, pvalue_table = TRUE)
ggsave("GSEA_top3.png", width = 9, height = 6)

ridgeplot(gse_h)   # opcional

write.csv(as.data.frame(gse_go), "GSEA_GO.csv", row.names = FALSE)
write.csv(as.data.frame(gse_h),  "GSEA_Hallmark.csv", row.names = FALSE)
```

## Arquivos de input e output
[Pasta no Google Drive](https://drive.google.com/drive/folders/1QDaOsWwIyT8FkGg632ImXi6VpDVscEBb?usp=drive_link)
