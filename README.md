
# CRISPR Pipeline - README

## Overview
 The pipeline includes downloading raw sequencing data, processing it, and identifying genetic variants.

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


(Step 3) During Quality Control check the following parameters:
1. Per base sequence quality
2. Overrepresented sequences
3. Adapter content


(Step 4) For Trimming
These are the Universal Adapters: 
_Illumina Universal Adapter_		=			AGATCGGAAGAG  <br />
_Illumina Small RNA 3' Adapter_		=		TGGAATTCTCGG  <br />
_Illumina Small RNA 5' Adapter_	=		GATCGTCGGACT  <br />
_Nextera Transposase Sequence_		=		CTGTCTCTTATA  <br />
_PolyA_			=		AAAAAAAAAAAA  <br />
_PolyG_		=		GGGGGGGGGGGG  <br />

Open newly created **adapter.fasta** file in Notepad and write:   <br />
```bash
>H1
AGATCGGAAGAG
```
