#!/bin/bash

#--------------------------------------------------------------------#
#==============DO NOT activate if launching with nohup===============#
#--------------------------------------------------------------------#
#spinner() {                                                         #
#  local pid=$1                                                      #
#  local message=$2                                                  #
#  local delay=0.1                                                   #
#  local spinstr='|/-\'                                              #
#  printf "%s " "$message"                                           #
#  while ps -p "$pid" > /dev/null; do                                #
#     local temp=${spinstr#?}                                        #
#      printf "[%c]  " "$spinstr"                                    #
#      spinstr=$temp${spinstr%"$temp"}                               #
#      sleep $delay                                                  #
#      printf "\b\b\b\b\b"                                           #
# done                                                               #
# printf "    \b\b\b\b"  # Clear spinner                             #
#}                                                                   #
#--------------------------------------------------------------------#

#----------------------------------------------------------------------------#
#==================Time function to measure execution time===================#
#----------------------------------------------------------------------------#
timed_exec() {
    local msg="$1"
    shift
    local start=$(date +%s)
    echo -n "START: $msg"
    "$@" &
    local cmd_pid=$!
    wait $cmd_pid
    local exit_code=$?
    local end=$(date +%s)
    local duration=$((end - start))
    echo -e "\nEND: $msg (Exit code: $exit_code) - Duration: ${duration}s"
    return $exit_code
}
#----------------------------------------------------------------------------#

#=============================== 1. CONFIGURATION ==============================#
#-----------------------------------------------------------------------------------------------------------------#
# Change the reference genome path if needed                                                                      #
# This is the reference genome used in the script                                                                 #
# It should be the same as the one used in the prepare_ref.sh script                                              #
# If you change it, make sure to update the prepare_ref.sh script accordingly                                     #
# Change the paths below if your tools are in different locations                                                 #
#-----------------------------------------------------------------------------------------------------------------#
# Reference genome                                                                                                #
REF="/home/usr/project/ref_genome/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"                                          #
# Path to FastQC binary                                                                                           #
FASTQC="/home/usr/project/Programs/FastQC/fastqc"                                                                 #
# Path to Adapter Removal                                                                                         #
FASTP="/home/usr/project/Programs/tools/fastp/fastp"                                                              #
# Path to BWA                                                                                                     #
BWA="/home/usr/project/Programs/tools/bwa-0.7.17/bwa"                                                             #
# Path to Samtools                                                                                                #
SAMTOOLS="/home/usr/project/Programs/tools/samtools-1.19.2/samtools"                                              #
# Path to GATK                                                                                                    #
GATK="/home/usr/project/Programs/tools/gatk-4.5.0.0/gatk"                                                         #
# Path to bcftools                                                                                                #
BCFTOOLS="/home/usr/project/Programs/tools/bcftools-1.19/bcftools"                                                #
# Path to Qualimap                                                                                                #
QUALIMAP="/home/usr/project/Programs/tools/qualimap_v2.3/qualimap"                                                #
#-----------------------------------------------------------------------------------------------------------------#

#=============================== 2. PATH SETUP =================================#
BASE_DIR="$(pwd)"
FASTQ_DIR="$BASE_DIR/raw_data"
OUTPUT_DIR="$BASE_DIR/output"
mkdir -p "$OUTPUT_DIR"
#============================================================================#

#============================== MAIN SCRIPT =================================#
shopt -s nullglob
FILES=("$FASTQ_DIR"/*_R1_*.fastq.gz)

if [ "${#FILES[@]}" -eq 0 ]; then
    echo "---------------------No FASTQ R1 found in $FASTQ_DIR-------------------"
    exit 1
fi

for R1 in "${FILES[@]}"; do
    R2="${R1/_R1_/_R2_}"
    NAME=$(basename "$R1" | sed 's/_R1_.*.fastq.gz//')
    SAMPLE_OUT="$OUTPUT_DIR/$NAME"
    FASTQC_OUT="$SAMPLE_OUT/fastqc_results"

    mkdir -p "$FASTQC_OUT"

    echo "======================3.FastQC on $NAME======================"
    timed_exec "FastQC on $NAME R1" "$FASTQC" -o "$FASTQC_OUT" -f fastq "$R1" --quiet
    timed_exec "FastQC on $NAME R2" "$FASTQC" -o "$FASTQC_OUT" -f fastq "$R2" --quiet

    echo "======================4.Fastp on $NAME======================"
    timed_exec "Fastp on $NAME" "$FASTP" \
       -i "$R1" -I "$R2" \
       -o "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" \
       -O "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" \
       --detect_adapter_for_pe \
       --trim_front1 20 --trim_front2 20 \
       --thread 4 --qualified_quality_phred 30 \
       --length_required 30 \
       --html "$SAMPLE_OUT/${NAME}_fastp.html" \
       --json "$SAMPLE_OUT/${NAME}_fastp.json"

    echo "==================5.Second FastQC analysis on $NAME=================="
    timed_exec "FastQC R1 trimmed $NAME" "$FASTQC" -o "$FASTQC_OUT" -f fastq "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" --quiet
    timed_exec "FastQC R2 trimmed $NAME" "$FASTQC" -o "$FASTQC_OUT" -f fastq "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" --quiet

    echo "==================6.BWA MEM on $NAME=================="
    "$BWA" mem -t 4 \
    -R "@RG\tID:${NAME}\tSM:${NAME}\tPL:illumina" \
    "$REF" \
    "$SAMPLE_OUT/${NAME}_R1_trimmed.fastq.gz" "$SAMPLE_OUT/${NAME}_R2_trimmed.fastq.gz" > "$SAMPLE_OUT/${NAME}.sam"

    echo "==================7.SAM to BAM conversion on $NAME=================="
    timed_exec "SAM to BAM $NAME" "$SAMTOOLS" view -Sb "$SAMPLE_OUT/${NAME}.sam" -o "$SAMPLE_OUT/${NAME}.bam"

    echo "==================8.Sorting BAM by name on $NAME=================="
    timed_exec "Sort BAM by name $NAME" "$SAMTOOLS" sort -n -m 3G "$SAMPLE_OUT/${NAME}.bam" -o "$SAMPLE_OUT/${NAME}_namesorted.bam"

    echo "==================8.Fixmate on $NAME=================="
    timed_exec "Fixmate $NAME" "$SAMTOOLS" fixmate -m "$SAMPLE_OUT/${NAME}_namesorted.bam" "$SAMPLE_OUT/${NAME}_fixmate.bam"

    echo "==================8.Sorting BAM by coordinate on $NAME=================="
    timed_exec "Sort BAM by coordinate $NAME" "$SAMTOOLS" sort -m 3G "$SAMPLE_OUT/${NAME}_fixmate.bam" -o "$SAMPLE_OUT/${NAME}_sorted.bam"

    echo "==================9.Mark duplicates on $NAME=================="
    timed_exec "Mark duplicates $NAME" "$SAMTOOLS" markdup -r "$SAMPLE_OUT/${NAME}_sorted.bam" "$SAMPLE_OUT/${NAME}_rmdup.bam"
    timed_exec "Index final BAM $NAME" "$SAMTOOLS" index "$SAMPLE_OUT/${NAME}_rmdup.bam"

    echo "==================10.Qualimap on $NAME=================="
    mkdir -p "$SAMPLE_OUT/qualimap"
    unset DISPLAY
    timed_exec "Qualimap $NAME" "$QUALIMAP" bamqc -bam "$SAMPLE_OUT/${NAME}_rmdup.bam" \
        -outdir "$SAMPLE_OUT/qualimap" -outformat PDF:HTML --java-mem-size=30G

    echo "==================11.HaplotypeCaller on $NAME=================="
    timed_exec "HaplotypeCaller $NAME" "$GATK" HaplotypeCaller \
        -R "$REF" \
        -I "$SAMPLE_OUT/${NAME}_rmdup.bam" \
        -O "$SAMPLE_OUT/${NAME}.vcf.gz" \
        -ERC GVCF

    if [ $? -eq 0 ]; then
       echo "VCF generated successfully, indexing..."
       timed_exec "Indexing VCF on $NAME" "$BCFTOOLS" index "$SAMPLE_OUT/${NAME}.vcf.gz"
       echo "Indexing done."
    else
       echo "Error generating VCF."
       exit 1
    fi

    echo "==================12.Cleaning temporary files=================="
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

    echo "----DONE $NAME----"
done
