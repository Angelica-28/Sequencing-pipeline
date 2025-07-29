#!/bin/bash

# === CONFIGURAZIONE ===
REF="Bos_taurus.ARS-UCD1.2.dna.toplevel.fa"  # <-- cambia se il tuo file ha nome diverso
GATK="/home/angelica/WGS/Programs/tools/gatk-4.5.0.0/gatk"
SAMTOOLS="/home/angelica/WGS/Programs/tools/samtools-1.19.2/samtools"
BWA="/home/angelica/WGS/Programs/tools/bwa-0.7.17/bwa"
# === CHECK FILE ===
if [ ! -f "$REF" ]; then
  echo "❌ Reference file $REF non trovato. Controlla il nome."
  exit 1
fi

echo "✅ Inizio indicizzazione della reference: $REF"

# 1. Samtools faidx
echo "→ Genero file .fai (samtools faidx)..."
$SAMTOOLS faidx "$REF"

# 2. Picard/GATK dictionary
echo "→ Genero file .dict (CreateSequenceDictionary)..."
$GATK CreateSequenceDictionary -R "$REF"

# 3. BWA index
echo "→ Genero index per BWA..."
$BWA index "$REF"

echo "✅ Tutti i file di indicizzazione generati con successo:"
echo "  - ${REF}.fai"
echo "  - ${REF%.fa}.dict"
echo "  - file .amb, .ann, .bwt, .pac, .sa (per BWA)"