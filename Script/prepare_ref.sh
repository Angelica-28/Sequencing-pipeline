#!/bin/bash

# === CONFIGURATION ===
REF="Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"  # <-- change if necessary
GATK="/home/usr/project/Programs/tools/gatk-4.5.0.0/gatk"
SAMTOOLS="/home/usr/project/Programs/tools/samtools-1.19.2/samtools"
BWA="/home/usr/project/Programs/tools/bwa-0.7.17/bwa"
# === CHECK FILE ===
if [ ! -f "$REF" ]; then
  echo "❌ Reference file $REF not found. Check name."
  exit 1
fi

echo "✅ Indexing: $REF"

# 1. Samtools faidx
echo "→ Generating file .fai (samtools faidx)..."
$SAMTOOLS faidx "$REF"

# 2. Picard/GATK dictionary
echo "→ Generating file .dict (CreateSequenceDictionary)..."
$GATK CreateSequenceDictionary -R "$REF"

# 3. BWA index
echo "→ Generating index per BWA..."
$BWA index "$REF"

echo "✅ All index file generated:"
echo "  - ${REF}.fai"
echo "  - ${REF%.fa}.dict"
echo "  - file .amb, .ann, .bwt, .pac, .sa (per BWA)"
