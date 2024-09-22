
# CRISPR Pipeline - README

## Overview

This repository provides a step-by-step guide for performing quality control, trimming, alignment, variant calling, and annotation of CRISPR data. The pipeline includes downloading raw sequencing data, processing it, and identifying genetic variants.

### Prerequisites

Ensure that the following tools are installed on your system:
- `wget`: Used to download raw sequencing data and reference genomes.
- `gunzip`: For decompressing `.gz` files.
- [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): For quality control of raw reads.
- [`fastp`](https://github.com/OpenGene/fastp): For trimming sequencing reads.
- [`BWA`](http://bio-bwa.sourceforge.net/): For aligning the sequencing reads to a reference genome.
- [`samtools`](http://www.htslib.org/): For sorting, removing duplicates, and converting file formats.
- [`GATK`](https://gatk.broadinstitute.org/hc/en-us): For variant calling.
- [`picard-tools`](http://broadinstitute.github.io/picard/): Required for file format conversion and preparing BAM files for GATK.
- [`SnpEff`](https://pcingola.github.io/SnpEff/): For variant annotation.
- [`VEP`](https://www.ensembl.org/info/docs/tools/vep/index.html): Used for annotating variants based on the genome.

## Pipeline Steps

### Step 1: Downloading Raw Sequencing Data
Use `wget` to download raw FASTQ files from the SRA repository:
```bash
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR217/020/SRR21763320/SRR21763320_1.fastq.gz
wget -c ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR217/020/SRR21763320/SRR21763320_2.fastq.gz
```

### Step 2: Data Extraction
Decompress the FASTQ files:
```bash
gunzip *.gz
```

### Step 3: Quality Control
Run quality control using `FastQC`:
```bash
fastqc *.fastq
```
Check the following parameters:
1. Per base sequence quality
2. Overrepresented sequences
3. Adapter content

### Step 4: Trimming
Trim the reads using `fastp`:
```bash
fastp -i sample.fastq -o trim_sample.fastq --adapter_fasta adapter.fasta
```

Run quality control on the trimmed reads:
```bash
fastqc trim_sample.fastq
```

### Step 5: Aligning Reads to the Reference Genome
Download the reference genome:
```bash
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
gunzip chr17.fa.gz
mv chr17.fa genome.fa
```

Index the reference genome:
```bash
bwa index -a bwtsw genome.fa
```

Align the reads to the reference genome:
```bash
bwa mem -t 2 genome.fa SRR21763320_1.fastq SRR21763320_2.fastq > bwa_SRR21763320.bam
```

### Step 6: Sorting and Conversion
Sort the BAM file using `samtools`:
```bash
samtools sort bwa_SRR21763320.bam > sorted_SRR21763320.bam
```

Convert the BAM file to SAM:
```bash
samtools view sorted_SRR21763320.bam > sorted_SRR21763320.sam
```

### Step 7: Removing Duplicates
Remove duplicate reads using `samtools`:
```bash
samtools rmdup -sS sorted_SRR21763320.bam rmdup_SRR21763320.bam
```

### Step 8: Variant Calling
Download and install GATK:
```bash
wget -c https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
```

Convert the reference genome to Picard-tools format:
```bash
picard-tools CreateSequenceDictionary R=genome.fa O=genome.dict
```

Prepare the BAM file for GATK:
```bash
picard-tools AddOrReplaceReadGroups I=rmdup_SRR21763320.bam O=picard_output.bam RGLB=lib1 RGPL=illumina RGPU=run RGSM=SRR21763320 SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
```

Call variants with GATK:
```bash
samtools faidx genome.fa
java -jar /path/to/gatk-package-4.3.0.0-local.jar HaplotypeCaller -R genome.fa -I picard_output.bam -O GATK_output.vcf
```

### Step 9: Variant Filtering
Filter variants using `SnpSift`:
```bash
cat GATK_output.vcf | java -jar /path/to/SnpSift.jar filter "(( QUAL>=30) & (DP>=10) & (MQ>=30))" > filter.vcf
```

### Step 10: Variant Annotation
Annotate variants using VEP:
Refer to the [VEP documentation](https://www.ensembl.org/info/docs/tools/vep/online/VEP_web_documentation.pdf) for more details.
```bash
java -jar snpEff.jar chr3 GATK_output.vcf > VEP_output.vcf
```

## License
This project is licensed under the MIT License.

## Acknowledgments
This pipeline was inspired by various bioinformatics tools and resources such as FastQC, BWA, GATK, and VEP.
