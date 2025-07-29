# 📦 FASTQ to PLINK-ready VCF Pipeline

This repository provides a full pipeline for processing raw sequencing data in FASTQ format and converting it into a VCF file that is compatible with PLINK. The pipeline is designed to work in a local environment without root privileges, using locally installed bioinformatics tools.

## 🧬 Overview

The pipeline includes the following steps:

1. **Quality control** with FastQC  
2. **Adapter trimming and quality filtering** using AdapterRemoval  
3. **Read alignment** to a reference genome using BWA  
4. **SAM to BAM conversion** using Samtools  
5. **Variant calling** with BCFtools mpileup  
7. **VCF preparation** for PLINK compatibility  

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
  - [`Bcftools`](https://samtools.github.io/bcftools/)
  - Java 8+ (required for GATK)
  - [`DeDup`](0.12.8) - Present in repository file
  - [`Qualimap`](http://qualimap.conesalab.org/) - Downloaded

## 📁 Folder structure

# The recommended directory organization is:
```plaintext
project/
│
├── raw_data/
│   ├── samples/# Input FASTQ files
│   ├── pipelinefromFastqtoBAM.sh
├── ref_genome/ # Reference genome FASTA and index files
│   ├── prepare_ref.sh
├── Programs
    ├── FastQC
    ├── tools/ # Locally installed bioinformatics tools
        ├── adapterremoval
        ├── bwa
        ├── bcftools
        ├── Dedup
        ├── gatk
        ├── qualimap
        ├── samtools
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

### 🔹 Bcftools
```bash
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -xjf bcftools-1.19.tar.bz2
cd bcftools-1.19
make
```
## 🔁 Pipeline steps explained

The script `pipelinefromFastqtoBAM.sh` automates the full preprocessing pipeline, from raw FASTQ files to final cleaned BAM files ready for variant calling and PLINK conversion. Here's a step-by-step explanation of what each section does:

1. **Setup and tool configuration**  
   The script sets up environment variables pointing to the reference genome and all required tools (FastQC, AdapterRemoval, BWA, Samtools, GATK, DeDup, Qualimap).
   ```bash
   # Reference genome
   REF="/home/usr/project/Ref_gen/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"

   # Tool paths (adjust these according to your system)
   FASTQC="/home/usr/project/Programs/FastQC/fastqc"
   ADAPTERREMOVAL="/home/usr/project/Programs/tools/adapterremoval/adapterremoval-2.3.4/build/AdapterRemoval"
   BWA="/home/usr/project/Programs/tools/bwa-0.7.17/bwa"
   SAMTOOLS="/home/usr/project/Programs/tools/samtools-1.19.2/samtools"
   GATK="/home/usr/project/Programs/tools/gatk-4.5.0.0/gatk"
   ```
> [!WARNING]
> These paths are specific to the local setup and must be edited before running the script. Make sure all tools are installed in accessible directories and update the variables accordingly.

3. **Folder navigation**  
   The script loops through all subfolders (one per sample), entering each to process the data found there.

4. **Initial FastQC**  
   Performs quality control on the raw FASTQ files using FastQC and stores reports in a dedicated folder.
   ```bash
   "$FASTQC" -o "$OUTPUT/$OUTPUT_DIR" -f fastq "$R1" --quiet
   "$FASTQC" -o "$OUTPUT/$OUTPUT_DIR" -f fastq "$R2" --quiet
   ```
6. **Adapter trimming**  
   Uses AdapterRemoval to clean the reads from adapter sequences, short reads (<30 bp), and low-quality bases.
   ```bash
   "$ADAPTERREMOVAL" \
            --file1 "../$R1" \
            --file2 "../$R2" \
            --basename "$NAME.adapt" \
            --minlength 30 \
            --minquality 30 \
            --trimns \
            --trimqualities \
            --collapse \
            --minalignmentlength 11 \
            --gzip
   ```
>[!NOTE]
> The `--collapse` option is used because we have paired-end reads (forward and reverse) and we want to merge overlapping pairs into a single read.  
7. **Second FastQC**  
   Re-runs FastQC on the cleaned (trimmed) reads to assess post-processing quality.
   ```bash
   "$FASTQC" -o "$OUTPUT/$OUTPUT_DIR" -f fastq "$OUTPUT/$NAME.adapt.collapsed.truncated.gz" --quiet
   ```
8. **Read alignment**  
   Aligns reads to the reference genome using `bwa aln`, and generates a `.sam` file with `bwa samse`.
   ```bash
   "$BWA" mem -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:illumina" \
   "$REF" \
   "$NAME.adapt.collapsed.fastq.gz" \
   > "$NAME.sam"
   ```
  >[!NOTE]
  > This step require time!      
10. **SAM to BAM conversion**  
   Converts the SAM file to BAM using `samtools view`.
   ```bash
  "$SAMTOOLS" view -Sb ${OUTPUT}/${NAME}.sam -o ${OUTPUT}/${NAME}.bam
   ```
11. **Sorting and indexing BAM**  
   Sorts and indexes the BAM file using `samtools sort` and `samtools index`.
   ```bash
   "$SAMTOOLS" sort -m 3G ${OUTPUT}/${NAME}.bam -o ${OUTPUT}/${NAME}_sort.bam
   "$SAMTOOLS" index ${OUTPUT}/${NAME}_sort.bam
   ```           
12. **Filtering by mapping quality (MAPQ ≥ 30)**  
   Filters out low-confidence alignments (MAPQ < 30), re-sorts, and re-indexes the BAM file.
   ```bash
   "$SAMTOOLS" view -q 30 -b ${OUTPUT}/${NAME}_sort.bam -o ${OUTPUT}/${NAME}_sort_mq.bam
   "$SAMTOOLS" sort -m 3G ${OUTPUT}/${NAME}_sort_mq.bam -o ${OUTPUT}/${NAME}_sort_mq_sorted.bam
   "$SAMTOOLS" index ${OUTPUT}/${NAME}_sort_mq_sorted.bam
   ```            
13. **Cleaning the BAM file**  
    Uses `gatk CleanSam` to fix potential formatting issues in the BAM file.
    ```bash
    gatk CleanSam -I ${OUTPUT}/${NAME}_sort_mq_sorted.bam -O ${OUTPUT}/${NAME}_sort_mq_sorted_clean.bam --VALIDATION_STRINGENCY SILENT 
    ```
>[!NOTE]
>you might have this error `/usr/bin/env: ‘python’: No such file or directory`
>To solve: `mkdir -p ~/bin`
>`ln -s "$(which python3)" ~/bin/python`
14. **Duplicate removal**  
    Removes PCR duplicates using `DeDup`. This improves variant calling quality.
    ```bash
    java -jar "$DEDUP" -i ${OUTPUT}/${NAME}_sort_mq_sorted_clean.bam -o ${OUTPUT} -m -u
    ```            
15. **Final BAM sorting and indexing**  
    The de-duplicated BAM is sorted and indexed again for downstream analysis.
    ```bash
    "$SAMTOOLS" sort -@ $CORES -m 3G ${OUTPUT}/${NAME}_sort_mq_sorted_clean_rmdup.bam -o ${OUTPUT}/${NAME}_rmdup_sort.bam
    "$SAMTOOLS" index ${OUTPUT}/${NAME}_rmdup_sort.bam
    ```            
16. **BAM quality summary**  
    Generates a `flagstat` summary with samtools to check alignment statistics.
    ```bash
    "$SAMTOOLS" flagstat ${OUTPUT}/${NAME}_rmdup_sort.bam > ${OUTPUT}/flagstat_${NAME}.txt
    ```            
17. **Qualimap report**  
    Produces a detailed HTML and PDF report with Qualimap for assessing alignment coverage and quality.
    ```bash
    mkdir -p ${OUTPUT}/Output/qualimap/
    unset DISPLAY 
    "$QUALIMAP" bamqc -bam ${OUTPUT}/${NAME}_rmdup_sort.bam -nt $CORES -outdir ${OUTPUT}/Output/qualimap/ -outformat PDF:HTML --java-mem-size=30G 
    ```            
18. **Variant calling with BCFtools**
    Generates a VCF file from the cleaned, sorted, and duplicate-removed BAM file using `bcftools mpileup` and `bcftools call`.
```bash
"$BCFTOOLS" mpileup -Ou -f "$REF" ${OUTPUT}/${NAME}_rmdup_sort.bam | \
"$BCFTOOLS" call -mv -Oz -o ${OUTPUT}/${NAME}.vcf.gz
"$BCFTOOLS" index ${OUTPUT}/${NAME}.vcf.gz
```
> [!WARNING]
> Make sure the reference genome is indexed using `samtools faidx`, `GATK CreateSequenceDictionary`, and `bwa index` (to generate `.fai`, `.dict`, and BWA index files).  
> A script is available in the repository as `prepare_ref.sh`.

## 🚀 Launching Script
Start your WGS pipeline with:

```bash
bash pipeline.sh
```
or 
```bash
nohup bash pipeline.sh > pipeline.log 2>&1 &
```
for background analysis and check status with ```tail pipeline.log```

## Output files 📁

The pipeline generates the following files for each processed sample inside the 📁 `${NAME}` folder:

- 📄 `${NAME}.sam` — raw alignment file produced by BWA MEM.
- 📄 `${NAME}.bam` — unsorted BAM file converted from SAM.
- 📂 `fastqc_results/` — directory containing FastQC reports.
- 🗂️ `Output/qualimap/` — directory containing Qualimap reports in PDF and HTML format.  
- 📝 `flagstat_${NAME}.txt` — mapping statistics report from `samtools flagstat`.
- 📄 `${NAME}.vcf.gz` — compressed VCF file containing the called variants.
- 📄 `${NAME}.vcf.gz.csi` — index file for the VCF.
