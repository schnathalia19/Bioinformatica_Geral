# Tarefa 1 — Do sequenciamento ao alinhamento

## Organismo

- **Nome científico:** *Mycoplasma genitalium*
- **Run:** `SRR39974648`
- **Dados:** reads paired-end

## Genoma de referência

- **Assembly:** `GCF_040556925.1_ASM4055692v1`
- **Arquivo FASTA:** `GCF_040556925.1_ASM4055692v1_genomic.fna`

### Download e conversão das reads

```bash
prefetch SRR39974648
fasterq-dump --split-files SRR39974648
```

### Indexação da referência

```bash
minimap2 -d M_genitalium.mmi GCF_040556925.1_ASM4055692v1_genomic.fna
```

### Alinhamento

```bash
minimap2 -ax sr M_genitalium.mmi SRR39974648_1.fastq SRR39974648_2.fastq > alinhamento_SRR39974648.sam
```

### Ordenação e indexação do BAM

```bash
samtools sort -o alinhamento_SRR39974648.bam alinhamento_SRR39974648.sam
samtools index alinhamento_SRR39974648.bam
```

## Arquivos de input e output

(https://drive.google.com/drive/folders/1pEP3cjxZ8WJUcEzd2Kppi2nfOzLeey7H?usp=sharing)
