# Tarefa 2 — Bulk RNA-seq

## Informações
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE348127
- SRA Run Selector: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1532893&o=acc_s%3Aa

## Download dos FASTQ
mkdir -p fastq
for r in SRR40776939 SRR40776940 SRR40776941 SRR40776942 SRR40776943 SRR40776944; do
  echo "Baixando $r..."
  fastq-dump -X 5000000 --split-3 --gzip --skip-technical --readids -O fastq "$r"
done

## Controle de qualidade
mkdir -p qc
fastqc fastq/*.fastq.gz -o qc
multiqc qc -o qc

## Índice do genoma e anotação
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar xzf grch38_genome.tar.gz
wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz
gunzip Homo_sapiens.GRCh38.112.gtf.gz

## Alinhamento + sort + index + flagstat
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

## Contagem
mkdir -p counts
featureCounts -p --countReadPairs -s 2 -T 8 \
  -a Homo_sapiens.GRCh38.112.gtf \
  -o counts/counts.txt bam/*.bam
cat counts/counts.txt.summary

## Arquivos de input e output

(https://drive.google.com/drive/folders/1QDaOsWwIyT8FkGg632ImXi6VpDVscEBb?usp=drive_link)
