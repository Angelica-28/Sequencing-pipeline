# 📦 FASTQ to PLINK-ready VCF Pipeline

This repository provides a full pipeline for processing raw sequencing data in FASTQ format and converting it into a VCF file that is compatible with PLINK. The pipeline is designed to work in a local environment without root privileges, using locally installed bioinformatics tools.

## 🧬 Overview

The pipeline includes the following steps:

1. **Quality control** with FastQC 
2. **Adapter trimming and quality filtering** using Fastp
3. **Read alignment** to a reference genome using BWA  
4. **SAM to BAM conversion** using Samtools
5. **Coverage control** with Qualimap
6. **Variant calling** with GATK HaplotypeCaller 
7. **VCF preparation** for PLINK compatibility  

This workflow is intended for whole-genome or exome sequencing data. It can be adapted for RNA-seq or targeted sequencing datasets.

## ⚙️ Requirements

To run the pipeline, you need:

- Bash shell (Unix-like environment)
- Local (non-root) installations of:
  - [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
  - [`Fastp`](http://opengene.org/fastp/fastp.1.0.1)
  - [`BWA`](http://bio-bwa.sourceforge.net/)
  - [`Samtools`](http://www.htslib.org/)
  - [`GATK`](https://gatk.broadinstitute.org/)
  - [`Bcftools`](https://samtools.github.io/bcftools/)
  - Java 8+ (required for GATK)
  - [`Qualimap`](http://qualimap.conesalab.org/) - Downloaded

## 📁 Folder structure

# The recommended directory organization is:
```plaintext
📂 project/
│
├── 📂 raw_data/
│   ├── 📂 samples/# Input FASTQ files
│   ├── 📝 seq_pipeline.sh
├── 📂 ref_genome/ # Reference genome FASTA and index files
│   ├── 📝 prepare_ref.sh
├──📂  Programs
    ├── 📂 FastQC
    ├── 📂 tools/ # Locally installed bioinformatics tools
        ├──📂  bwa
        ├── 📂 bcftools
        ├── 📂 gatk
        ├── 📂 qualimap
        ├── 📂 samtools
        ├── 📂 fastp
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

### 🔹 Fastp

```bash
mkdir -p tools/fastp
cd tools/fastqc
wget http://opengene.org/fastp/fastp.1.0.1
chmod a+x fastp.1.0.1
mv fastp.1.0.1 fastp
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
mkdir -p tools/bcftools
cd tools/bcftools
wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
tar -xjf bcftools-1.19.tar.bz2
cd bcftools-1.19
make
```

### 🔹 Picard
```bash
mkdir -p tools/picard
cd tools/picard
wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar
```
## 🔁 Pipeline steps explained

The script `seq_pipeline.sh` automates the full preprocessing pipeline, from raw FASTQ files to final cleaned BAM files ready for variant calling and PLINK conversion. Here's a step-by-step explanation of what each section does:

1. **Setup and tool configuration**  
   The script sets up environment variables pointing to the reference genome and all required tools (FastQC, AdapterRemoval, BWA, Samtools, GATK, DeDup, Qualimap).
   ```bash
   # Reference genome
   REF="/home/usr/project/Ref_gen/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"
   FASTQC="/home/usr/WGS/Programs/FastQC/fastqc"                                                                                                                                                     
   FASTP="/home/usr/WGS/Programs/tools/fastp/fastp"                                                                                                                                                              
   BWA="/home/usr/WGS/Programs/tools/bwa-0.7.17/bwa"                                                                                                                                                          
   SAMTOOLS="/home/usr/WGS/Programs/tools/samtools-1.19.2/samtools"                                                                                                                                                  
   GATK="/home/usr/WGS/Programs/tools/gatk-4.5.0.0/gatk"                                                                                                                                                       
   BCFTOOLS="/home/usr/WGS/Programs/tools/bcftools-1.19/bcftools"                                                                                                                                                                                    
   QUALIMAP="/home/usr/WGS/Programs/tools/qualimap_v2.3/qualimap"                                                                                             
   PICARD="/home/angelica/WGS/Programs/tools/picard/picard.jar"  
   ```
> [!WARNING]
> These paths are specific to the local setup and must be edited before running the script. Make sure all tools are installed in accessible directories and update the variables accordingly.

2. **Folder navigation**  
   The script loops through all subfolders (one per sample), entering each to process the data found there.

3. **Initial FastQC**  
   Performs quality control on the raw FASTQ files using FastQC and stores reports in a dedicated folder.
   ```bash
   "$FASTQC" -o "$OUTPUT/$OUTPUT_DIR" -f fastq "$R1" --quiet
   "$FASTQC" -o "$OUTPUT/$OUTPUT_DIR" -f fastq "$R2" --quiet
   ```
4. **Fastp paired-end read trimming and quality filtering**  
   Fastp performs automatic detection and removal of adapter sequences in paired-end reads. In this configuration, the first 20 bases of each read in R1 and R2 are trimmed to mitigate potential biases or low-quality bases introduced during the initial sequencing cycles. Bases are classified as qualified only if their Phred score is ≥ 30 (corresponding to an estimated base-calling accuracy of 99.9%), which defines the threshold applied in subsequent quality filtering steps. Following trimming and quality control, reads shorter than 30 bp are discarded from the output.
   ```bash
       "$FASTP" \
       -i "$R1" -I "$R2" \
       -o "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" \
       -O "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" \
       --detect_adapter_for_pe \
       --trim_front1 20 --trim_front2 20 \
       --thread 4 --qualified_quality_phred 30 \
       --length_required 30 \
       --html "$SAMPLE_OUT/${NAME}_fastp.html" \
       --json "$SAMPLE_OUT/${NAME}_fastp.json"
   ``` 
5. **Second FastQC**  
   Re-runs FastQC on the cleaned (trimmed) reads to assess post-processing quality.
   ```bash
   "$FASTQC" -o "$FASTQC_OUT" -f fastq "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" --quiet
   "$FASTQC" -o "$FASTQC_OUT" -f fastq "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" --quiet
   ```
6. **Read alignment**  
   Aligns reads to the reference genome using `bwa mem`.
   ```bash
   "$BWA" mem -t 4 \
    -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:illumina" \
    "$REF" \
    "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" > "$SAMPLE_OUT/${NAME}.sam"
   ```

> [!WARNING]
> Make sure the reference genome is indexed using `samtools faidx`, `GATK CreateSequenceDictionary`, and `bwa index` (to generate `.fai`, `.dict`, and BWA index files).  
> A script is available in the repository as `prepare_ref.sh`.

>[!NOTE]
> This step require time!      

7. **SAM to BAM conversion**  
   Converts the SAM file to BAM using `samtools view`.
   ```bash
   "$SAMTOOLS" view -Sb ${OUTPUT}/${NAME}.sam -o ${OUTPUT}/${NAME}.bam
   ```
8. **Sorting and indexing BAM**  
   Samtools is used to sort and fix paired-end alignment information. The BAM file is first name-sorted, then mate information is fixed, and finally the BAM is coordinate-sorted to produce the final `${NAME}_sorted.bam` ready for downstream analysis.
   ```bash
   "$SAMTOOLS" sort -n -m 3G "$SAMPLE_OUT/${NAME}.bam" -o "$SAMPLE_OUT/${NAME}_namesorted.bam"
   "$SAMTOOLS" fixmate -m "$SAMPLE_OUT/${NAME}_namesorted.bam" "$SAMPLE_OUT/${NAME}_fixmate.bam"
   "$SAMTOOLS" sort -m 3G "$SAMPLE_OUT/${NAME}_fixmate.bam" -o "$SAMPLE_OUT/${NAME}_sorted.bam"
   ```
>[!NOTE]
> Coordinate-sorted BAM is required for duplicate marking.

9. **Duplicate removal**  
    Picard is used to mark and remove duplicates, then index the resulting BAM file for downstream analysis.
   ```bash
    java -jar "$PICARD" MarkDuplicates -I "$SAMPLE_OUT/${NAME}_sorted.bam" -O "$SAMPLE_OUT/${NAME}_rmdup.bam" -M "$SAMPLE_OUT/${NAME}_markdup_metrics.txt" -REMOVE_DUPLICATES true -CREATE_INDEX true
   ```
10. **Filter for mapping quality**

    Filter out any alignments with a mapping quality below 25, saving a BAM file that contains only high-quality, uniquely aligned reads. The file is then indexed for fast access and sorted to ensure proper coordinate order, preparin
    it for downstream variant calling and qualimap.
    ```bash
    "$SAMTOOLS" view -b -q 25 "$SAMPLE_OUT/${NAME}_rmdup.bam" > "$SAMPLE_OUT/${NAME}_rmdup_mapq25.bam"
    "$SAMTOOLS" sort -m 3G "$SAMPLE_OUT/${NAME}_rmdup_mapq25.bam" -o "$SAMPLE_OUT/${NAME}_rmdup_mapq25_sorted.bam"
    "$SAMTOOLS" index "$SAMPLE_OUT/${NAME}_rmdup_mapq25_sorted.bam"
    ```
12. **Qualimap report**  
    Produces a detailed HTML and PDF report with Qualimap for assessing alignment coverage and quality.
    ```bash
    mkdir -p "$SAMPLE_OUT/qualimap"
    unset DISPLAY
    timed_exec "Qualimap $NAME" "$QUALIMAP" bamqc -bam "$SAMPLE_OUT/${NAME}_rmdup_mapq25_sorted.bam" \
        -outdir "$SAMPLE_OUT/qualimap" -outformat PDF:HTML --java-mem-size=30G
    ```            
13. **Variant calling with GATK**

    Generates a VCF file from the cleaned, sorted, and duplicate-removed BAM file using `GATK HaplotypeCallet`.
    ```bash
    "$GATK" HaplotypeCaller \
        -R "$REF" \
        -I "$SAMPLE_OUT/${NAME}_rmdup.bam" \
        -O "$SAMPLE_OUT/${NAME}.vcf.gz" \
        -ERC GVCF
    ```
15. **Temp file cleaning**
    ```bash
    rm -f "$SAMPLE_OUT/${NAME}.sam" \
          "$SAMPLE_OUT/${NAME}.bam" \
          "$SAMPLE_OUT/${NAME}.bam.bai" \
          "$SAMPLE_OUT/${NAME}_fixmate.bam" \
          "$SAMPLE_OUT/${NAME}_fixmate.bam.bai" \
          "$SAMPLE_OUT/${NAME}_namesorted.bam" \
          "$SAMPLE_OUT/${NAME}_namesorted.bam.bai" \
          "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" \
          "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" \
          "$SAMPLE_OUT/${NAME}_sorted.bam" \
          "$SAMPLE_OUT/${NAME}_sorted.bam.bai" \
          "$R1" \
          "$R2"
     ```
> [!WARNING]
> original FASTQ files are removed in this step, since a copy is used; the originals are stored locally.

## 🚀 Launching Script
Start your WGS pipeline with:

```bash
bash seq_pipeline.sh
```
or 
```bash
nohup bash seq_pipeline.sh > pipeline.log 2>&1 &
```
for background analysis and check status with ```tail pipeline.log```

## Output files 📁

The pipeline generates the following files for each processed sample inside the 📁 `${NAME}` folder inside the 📁 `output` folder:

- 📄 `${NAME}.sam` — raw alignment file produced by BWA MEM.
- 📄 `${NAME}.bam` — unsorted BAM file converted from SAM.
- 📂 `fastqc_results/` — directory containing FastQC reports.
- 🗂️ `Output/qualimap/` — directory containing Qualimap reports in PDF and HTML format.  
- 📝 `flagstat_${NAME}.txt` — mapping statistics report from `samtools flagstat`.
- 📄 `${NAME}.vcf.gz` — compressed VCF file containing the called variants.
- 📄 `${NAME}.vcf.gz.csi` — index file for the VCF.
