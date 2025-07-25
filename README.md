# Sequencing-pipeline
This repository provides a full pipeline for processing raw sequencing data in FASTQ format and converting it into a VCF file that is compatible with PLINK. The pipeline is designed to work in a local environment without root privileges, using locally installed bioinformatics tools.
# 📦 FASTQ to PLINK-ready VCF Pipeline

This repository provides a full pipeline for processing raw sequencing data in FASTQ format and converting it into a VCF file that is compatible with PLINK. The pipeline is designed to work in a local environment without root privileges, using locally installed bioinformatics tools.

## 🧬 Overview

The pipeline includes the following steps:

1. **Quality control** with FastQC  
2. **Adapter trimming and quality filtering** using AdapterRemoval  
3. **Read alignment** to a reference genome using BWA  
4. **SAM to BAM conversion** using Samtools  
5. **Variant calling** with GATK HaplotypeCaller  
6. **VCF preparation** for PLINK compatibility  

This workflow is intended for whole-genome or exome sequencing data. It can be adapted for RNA-seq or targeted sequencing datasets.

## ⚙️ Requirements

To run the pipeline, you need:

- Bash shell (Unix-like environment)
- Local (non-root) installations of:
  - [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [`AdapterRemoval`](https://adapterremoval.readthedocs.io/en/stable/)
  - [`BWA`](http://bio-bwa.sourceforge.net/)
  - [`Samtools`](http://www.htslib.org/)
  - [`GATK`](https://gatk.broadinstitute.org/)
  - Java 8+ (required for GATK)
  - [`DeDup`](0.12.8) - Present in repository file
  - [`Qualimap`](v2.3) - Present in repository file

## 📁 Folder structure

# The recommended directory organization is:
```plaintext
project/
│
├── raw_data/# Input FASTQ files
  ├──Script.sh
├── ref_genome/ # Reference genome FASTA and index files
├── tools/ # Locally installed bioinformatics tools
  ├──FastQC
  ├──adapterremoval
  ├──bwa
  ├──Dedup
  ├──gatk
  ├──qualimap
  ├──samtools
```
## 🛠️ Installation

Since the pipeline is designed for environments without `sudo` access, all tools must be installed locally in the `tools/` directory.

You can download and install each tool with the following instructions:

### 🔹 FastQC

```bash
mkdir -p tools/fastqc
cd tools/fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod +x FastQC/fastqc
```
### 🔹 AdapterRemoval
```bash
mkdir -p tools/adapterremoval
cd tools/adapterremoval
wget https://github.com/MikkelSchubert/adapterremoval/archive/refs/tags/v2.3.4.tar.gz
tar -xvzf v2.3.4.tar.gz
cd adapterremoval-2.3.4
mkdir build && cd build
cmake ..
make
```

### 🔹 BWA
```bash
mkdir -p tools/bwa
cd tools/bwa
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
tar -xvjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
```

### 🔹 Samtools
```bash
mkdir -p tools/samtools
cd tools/samtools
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xvjf samtools-1.18.tar.bz2
cd samtools-1.18
./configure --prefix=$PWD/install
make
make install
```

### 🔹 GATK
```bash
mkdir -p tools/gatk
cd tools/gatk
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
```
## 🔁 Pipeline steps explained

The script `fastq_to_vcf.sh` automates the full preprocessing pipeline, from raw FASTQ files to final cleaned BAM files ready for variant calling and PLINK conversion. Here's a step-by-step explanation of what each section does:

1. **Setup and tool configuration**  
   The script sets up environment variables pointing to the reference genome and all required tools (FastQC, AdapterRemoval, BWA, Samtools, GATK, DeDup, Qualimap).
   ```bash
   # Reference genome
   REF="/home/angelica/WGS/Ref_gen/GCF_002263795.3_ARS-UCD2.0_cds_from_genomic.fna"

   # Tool paths (adjust these according to your system)
   FASTQC="/home/angelica/WGS/Programs/FastQC/fastqc"
   ADAPTERREMOVAL="/home/angelica/WGS/Programs/tools/adapterremoval/adapterremoval-2.3.4/build/AdapterRemoval"
   BWA="/home/angelica/WGS/Programs/tools/bwa-0.7.17/bwa"
   SAMTOOLS="/home/angelica/WGS/Programs/tools/samtools-1.19.2/samtools"
   GATK="/home/angelica/WGS/Programs/tools/gatk-4.5.0.0/gatk"
```

3. **Folder navigation**  
   The script loops through all subfolders (one per sample), entering each to process the data found there.

4. **Initial FastQC**  
   Performs quality control on the raw FASTQ files using FastQC and stores reports in a dedicated folder.

5. **Adapter trimming**  
   Uses AdapterRemoval to clean the reads from adapter sequences, short reads (<30 bp), and low-quality bases.

6. **Second FastQC**  
   Re-runs FastQC on the cleaned (trimmed) reads to assess post-processing quality.

7. **Read alignment**  
   Aligns reads to the reference genome using `bwa aln`, and generates a `.sam` file with `bwa samse`.

8. **SAM to BAM conversion**  
   Converts the SAM file to BAM using `samtools view`.

9. **Sorting and indexing BAM**  
   Sorts and indexes the BAM file using `samtools sort` and `samtools index`.

10. **Filtering by mapping quality (MAPQ ≥ 30)**  
   Filters out low-confidence alignments (MAPQ < 30), re-sorts, and re-indexes the BAM file.

11. **Cleaning the BAM file**  
    Uses `gatk CleanSam` to fix potential formatting issues in the BAM file.

12. **Duplicate removal**  
    Removes PCR duplicates using `DeDup`. This improves variant calling quality.

13. **Final BAM sorting and indexing**  
    The de-duplicated BAM is sorted and indexed again for downstream analysis.

14. **BAM quality summary**  
    Generates a `flagstat` summary with samtools to check alignment statistics.

15. **Qualimap report**  
    Produces a detailed HTML and PDF report with Qualimap for assessing alignment coverage and quality.

---

Let me know if vuoi anche aggiungere la sezione successiva, cioè il **variant calling** con GATK e la conversione del VCF in formato PLINK.
