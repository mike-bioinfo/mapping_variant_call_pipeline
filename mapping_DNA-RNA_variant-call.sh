#!/bin/bash

# ================================
# Clonal Pipeline for WES and RNA-seq data mapping and variant calling
# Author: Michael Torres
# Description:
#  
#   Requires: BWA, Samtools, GATK, featureCounts, R (DESeq2)
#   Note: This script requires a "SAMPLE_MAP" for renaming samples and 
#   correlate conditions during RNA mapping.
# ================================

set -euo pipefail

start_time=$(date +%s)

# # # Define colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color (reset)


# # # === CONFIGURATION ===
THREADS=18
REFERENCE="/path/to/databases/human-genome-ref/GRCh38.p14_genomic.primary.fna"
GTF="/path/to/databases/human-genome-ref/GRCh38.p14_genomic_primary.modified.gtf"
REF_GENOME_DIR="/path/to/databases/human-genome-ref/"
PON="/path/to/databases/ponDB/1000g_pon.hg38.modified.vcf.gz"
GNOMAD="/path/to/databases/gnomadDB/af-only-gnomad.hg38.vcf.modified.gz"
SAMPLE_MAP="/path/to/working-directory/sample_map.tsv"
FASTQ_DIR="/path/to/working-directory/data/fastq"
OUTDIR="/path/to/working-directory/results"
DBSNP="/path/to/databases/dbsnpDB/00-All.vcf.gz"

mkdir -p "$OUTDIR"/{bam,vcf,vcf/filtered,pyclone,rna_expression,plots,tmp}

# # # === READ SAMPLE MAP ===
declare -A SAMPLE_TYPE
declare -A WES_R1
declare -A WES_R2
declare -A RNA_REPS
declare -A PAIR_ID

SAMPLES=()

# # # Skip header line explicitly and read with proper column order
while IFS=$'\t' read -r sample_type sample_name pair_id wes_r1 wes_r2 rna_reps; do
    # Skip commented/empty lines and header
    [[ "$sample_type" =~ ^#.*$ ]] && continue
    [[ -z "$sample_name" ]] && continue
    
    # # Debug: Show raw input
    # echo "DEBUG_READ: type:$sample_type | name:$sample_name | pair:$pair_id | r1:$wes_r1 | r2:$wes_r2 | reps:$rna_reps"
    
    # # Assign to arrays - note corrected variable order
    SAMPLE_TYPE["$sample_name"]="$sample_type"
    WES_R1["$sample_name"]="$FASTQ_DIR/$wes_r1"
    WES_R2["$sample_name"]="$FASTQ_DIR/$wes_r2"
    RNA_REPS["$sample_name"]="$rna_reps"
    PAIR_ID["$sample_name"]="$pair_id"
    SAMPLES+=("$sample_name")
done < <(tail -n +2 "$SAMPLE_MAP")  # Skip first line (header)

# # === DEBUG: PRINT SAMPLE MAP VARIABLES ===
echo -e "\n${YELLOW}=== DEBUG: Variables from sample_map.tsv ===${NC}"

# Print header
printf "%-15s | %-10s | %-10s | %-30s | %-30s | %-50s\n" \
       "Sample" "Type" "Pair_ID" "WES_R1" "WES_R2" "RNA_Replicates"
echo "------------------------------------------------------------------------------------------------------------------------"

# Print each sample's information
for sample in "${SAMPLES[@]}"; do
    printf "%-15s | %-10s | %-10s | %-30s | %-30s | %-50s\n" \
           "$sample" \
           "${SAMPLE_TYPE[$sample]}" \
           "${PAIR_ID[$sample]}" \
           "${WES_R1[$sample]}" \
           "${WES_R2[$sample]}" \
           "${RNA_REPS[$sample]}"
done

# # Additional checks
echo -e "\n${BLUE}=== Array Sizes ===${NC}"
echo "Total samples: ${#SAMPLES[@]}"
echo "SAMPLE_TYPE entries: ${#SAMPLE_TYPE[@]}"
echo "PAIR_ID entries: ${#PAIR_ID[@]}"
echo "WES_R1 entries: ${#WES_R1[@]}"
echo "WES_R2 entries: ${#WES_R2[@]}"
echo "RNA_REPS entries: ${#RNA_REPS[@]}"

# # Verify specific samples
echo -e "\n${BLUE}=== Sample Verification ===${NC}"
check_samples=("HCA46_T0" "HCA46_CTX-R2" "HT29_T0" "HT29_DAB-CTX")
for s in "${check_samples[@]}"; do
    if [[ -z "${SAMPLE_TYPE[$s]}" ]]; then
        echo -e "${RED}ERROR: Sample $s not found in arrays!${NC}"
    else
        echo -e "${GREEN}Found sample $s:${NC}"
        echo "  Type: ${SAMPLE_TYPE[$s]}"
        echo "  Pair: ${PAIR_ID[$s]}"
        echo "  WES_R1: ${WES_R1[$s]}"
        echo "  WES_R2: ${WES_R2[$s]}"
        echo "  RNA_Reps: ${RNA_REPS[$s]}"
    fi
done

echo -e "\n${GREEN}=== Sample Map Loading Complete ===${NC}"

# # # === ALIGN WES DATA ===
echo -e "${GREEN}=== [Step 1] Aligning WES data... ===${NC}"

for sample in "${SAMPLES[@]}"; do
    echo " - Processing $sample"

    # Extract metadata
    sample_type="${SAMPLE_TYPE[$sample]}"
    pair_id="${PAIR_ID[$sample]}"
    RG="@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tST:${sample_type}"
    echo "   Read Group: $RG"
    
    bwa mem -t $THREADS -R "$RG" "$REFERENCE" "${WES_R1[$sample]}" "${WES_R2[$sample]}" | \
        samtools view -@ $THREADS -Sb - | \
        samtools sort -@ $THREADS -o "$OUTDIR/bam/${sample}.sorted.bam"
    samtools index "$OUTDIR/bam/${sample}.sorted.bam"

    gatk MarkDuplicates \
        -I "$OUTDIR/bam/${sample}.sorted.bam" \
        -O "$OUTDIR/bam/${sample}.dedup.bam" \
        -M "$OUTDIR/bam/${sample}.metrics.txt"
    
    gatk BaseRecalibrator -R "$REFERENCE" -I "$OUTDIR/bam/${sample}.dedup.bam" \
    --known-sites "$GNOMAD" --known-sites "$DBSNP" -O "$OUTDIR/bam/${sample}.recal.table"

    gatk ApplyBQSR -R "$REFERENCE" -I "$OUTDIR/bam/${sample}.dedup.bam" --bqsr-recal-file "$OUTDIR/bam/${sample}.recal.table" \
    -O "$OUTDIR/bam/${sample}.bqsr.bam"

    samtools index "$OUTDIR/bam/${sample}.bqsr.bam"
done

# # # === VARIANT CALLING ===
echo -e "${GREEN}=== [Step 2] Somatic mutation calling with Mutect2... ===${NC}"

for tumor in "${SAMPLES[@]}"; do
    if [[ "${SAMPLE_TYPE[$tumor]}" != "resist" ]]; then continue; fi
    for normal in "${SAMPLES[@]}"; do
        if [[ "${SAMPLE_TYPE[$normal]}" == "parent" && "${PAIR_ID[$tumor]}" == "${PAIR_ID[$normal]}" ]]; then
            pair_name="${tumor}_vs_${normal}"
            echo " - Calling variants for $pair_name"

            gatk Mutect2 \
                -R "$REFERENCE" \
                -I "$OUTDIR/bam/${tumor}.bqsr.bam" \
                -I "$OUTDIR/bam/${normal}.bqsr.bam" \
                --germline-resource "$GNOMAD" \
                --panel-of-normals "$PON" \
                -tumor "$tumor" \
                -normal "$normal" \
                -O "$OUTDIR/vcf/${pair_name}.vcf.gz"

            gatk FilterMutectCalls \
                -V "$OUTDIR/vcf/${pair_name}.vcf.gz" \
                -R "$REFERENCE" \
                -O "$OUTDIR/vcf/filtered/${pair_name}.filtered.vcf.gz"
        fi
    done
done

echo -e "${GREEN}=== [Step 2.1] Somatic mutation calling with Mutect2 for each tumor sample... ===${NC}"

for tumor in "${SAMPLES[@]}"; do
    if [[ "${SAMPLE_TYPE[$tumor]}" != "resist" ]]; then continue; fi
    for normal in "${SAMPLES[@]}"; do
        if [[ "${SAMPLE_TYPE[$normal]}" == "parent" && "${PAIR_ID[$tumor]}" == "${PAIR_ID[$normal]}" ]]; then
            pair_name="${tumor}_vs_${normal}"
            echo " - Calling variants for $tumor"

            gatk Mutect2 \
                -R "$REFERENCE" \
                -I "$OUTDIR/bam/${tumor}.bqsr.bam" \
                --germline-resource "$GNOMAD" \
                --panel-of-normals "$PON" \
                -O "$OUTDIR/vcf/${tumor}.vcf.gz"

            gatk FilterMutectCalls \
                -V "$OUTDIR/vcf/${tumor}.vcf.gz" \
                -R "$REFERENCE" \
                -O "$OUTDIR/vcf/filtered/${tumor}.filtered.vcf.gz"
        fi
    done
done

echo -e "${GREEN}=== [Step 2.2] Somatic mutation calling with Mutect2 for each normal sample... ===${NC}"

for tumor in "${SAMPLES[@]}"; do
    if [[ "${SAMPLE_TYPE[$tumor]}" != "resist" ]]; then continue; fi
    for normal in "${SAMPLES[@]}"; do
        if [[ "${SAMPLE_TYPE[$normal]}" == "parent" && "${PAIR_ID[$tumor]}" == "${PAIR_ID[$normal]}" ]]; then
            pair_name="${tumor}_vs_${normal}"
            echo " - Calling variants for $normal"

            gatk Mutect2 \
                -R "$REFERENCE" \
                -I "$OUTDIR/bam/${normal}.bqsr.bam" \
                --germline-resource "$GNOMAD" \
                --panel-of-normals "$PON" \
                -O "$OUTDIR/vcf/${normal}.vcf.gz"

            gatk FilterMutectCalls \
                -V "$OUTDIR/vcf/${normal}.vcf.gz" \
                -R "$REFERENCE" \
                -O "$OUTDIR/vcf/filtered/${normal}.filtered.vcf.gz"
        fi
    done
done

# # # === RNA-SEQ ALIGNMENT ===
echo -e "${GREEN}=== [Step 3] Aligning RNA-seq reads... ===${NC}"

# Create directory for STAR outputs
mkdir -p "$OUTDIR/rna_alignment"

for sample in "${SAMPLES[@]}"; do
    IFS=',' read -ra reps <<< "${RNA_REPS[$sample]}"
    for i in "${!reps[@]}"; do
        fq="$FASTQ_DIR/${reps[$i]}"
        out_prefix="$OUTDIR/rna_alignment/${sample}_rep$((i+1))_"
        
        echo -e "${YELLOW}=== - Aligning ${sample} replicate $((i+1)): ${reps[$i]} ===${NC}"
        
        STAR --runThreadN $THREADS \
             --genomeDir "$REF_GENOME_DIR" \
             --readFilesIn "$fq" \
             --readFilesCommand zcat \
             --outFileNamePrefix "$out_prefix" \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile "$GTF" \
             --outBAMcompression 6 \
             --outFilterType BySJout \
             --outFilterMultimapNmax 20 \
             --alignSJoverhangMin 8 \
             --alignSJDBoverhangMin 1 \
             --outFilterMismatchNmax 999 \
             --outFilterMismatchNoverLmax 0.04 \
             --alignIntronMin 20 \
             --alignIntronMax 1000000 \
             --alignMatesGapMax 1000000
        
        # Rename, merge and index the sorted BAM file from RNA-seq
        mv "${out_prefix}Aligned.sortedByCoord.out.bam" "$OUTDIR/bam/${sample}_rep$((i+1)).sorted.bam"
        samtools index "$OUTDIR/bam/${sample}_rep$((i+1)).sorted.bam"
        echo -e "${GREEN}=== Completed alignment for ${sample} replicate $((i+1)) ===${NC}"
    done
done

# # # === MERGING REPLICATES ===

for sample in "${SAMPLES[@]}"; do
    IFS=',' read -ra reps <<< "${RNA_REPS[$sample]}"
    bam_files=()  # Array to store BAM files for merging

    echo -e "${YELLOW}=== Processing ${sample} (${#reps[@]} replicates) ===${NC}"

    # Collect all BAM files for the sample
    for i in "${!reps[@]}"; do
        bam_files+=("$OUTDIR/bam/${sample}_rep$((i+1)).sorted.bam")
    done

    # Merge all replicates at once (if >1 replicate)
    if [[ ${#bam_files[@]} -gt 1 ]]; then
        echo -e "${YELLOW}Merging ${#bam_files[@]} replicates...${NC}"
        samtools merge -f "$OUTDIR/rna_alignment/${sample}_rna_merge.bam" "${bam_files[@]}"
    else
        echo -e "${YELLOW}Only 1 replicate; copying instead of merging...${NC}"
        cp "${bam_files[0]}" "$OUTDIR/rna_alignment/${sample}_rna_merge.bam"
    fi

    # Sort and index the merged BAM (once per sample)
    samtools sort -@ $THREADS "$OUTDIR/rna_alignment/${sample}_rna_merge.bam" -o "$OUTDIR/rna_alignment/${sample}_rna_sorted.bam"
    samtools index "$OUTDIR/rna_alignment/${sample}_rna_sorted.bam"

    echo -e "${GREEN}=== Processed ${sample} located at "$OUTDIR/rna_alignment/${sample}_rna_sorted.bam" ===${NC}"

done

# # #=== EXPRESSION QUANTIFICATION ===

echo -e "${GREEN}=== [Step 4] Running featureCounts... ===${NC}"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR/rna_expression"

# Initialize array for BAM files
BAM_FILES=()

# Collect all RNA-seq BAM files
for sample in "${SAMPLES[@]}"; do
    # Find all replicate BAMs for this sample
    replicates=("$OUTDIR/bam/${sample}_rep"*.sorted.bam)
    
    if [ ${#replicates[@]} -eq 0 ]; then
        echo -e "${RED}ERROR: No BAM files found for sample $sample${NC}"
        echo "Expected files like: $OUTDIR/bam/${sample}_rep1.sorted.bam"
        exit 1
    fi

    for bam in "${replicates[@]}"; do
        if [[ ! -f "$bam" ]]; then
            echo -e "${RED}ERROR: BAM file not found: $bam${NC}"
            exit 1
        fi
        BAM_FILES+=("$bam")
        echo " - Found BAM: $bam"
    done
done

# Check if any BAM files were found
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo -e "${RED}ERROR: No BAM files found for feature counting${NC}"
    exit 1
fi

echo -e "\n${BLUE}Running featureCounts on ${#BAM_FILES[@]} BAM files...${NC}"

# Get unique pair IDs
declare -A unique_pairs
for sample in "${SAMPLES[@]}"; do
    unique_pairs["${PAIR_ID[$sample]}"]=1
done

# Process each pair separately
for pair_id in "${!unique_pairs[@]}"; do
    echo -e "\n${BLUE}Processing pair: $pair_id${NC}"
    
    # Initialize arrays for parent and resistant BAMs
    PARENT_BAMS=()
    RESIST_BAMS=()

    # Collect BAM files for this pair
    for sample in "${SAMPLES[@]}"; do
        if [[ "${PAIR_ID[$sample]}" == "$pair_id" ]]; then
            # Find all replicate BAMs for this sample
            replicates=("$OUTDIR/bam/${sample}_rep"*.sorted.bam)
            
            if [[ "${SAMPLE_TYPE[$sample]}" == "parent" ]]; then
                PARENT_BAMS+=("${replicates[@]}")
            elif [[ "${SAMPLE_TYPE[$sample]}" == "resist" ]]; then
                RESIST_BAMS+=("${replicates[@]}")
            fi
        fi
    done

    # Run featureCounts for parent samples of this pair
    if [ ${#PARENT_BAMS[@]} -gt 0 ]; then
        echo -e "${YELLOW}Counting parent samples (${#PARENT_BAMS[@]} BAMs)${NC}"
        featureCounts -T $THREADS \
            -a "$GTF" \
            -o "$OUTDIR/rna_expression/${pair_id}_parent_counts.txt" \
            -g gene_id \
            -t exon \
            --countReadPairs \
            -s 1 \
            "${PARENT_BAMS[@]}"
    else
        echo -e "${RED}WARNING: No parent BAMs found for pair $pair_id${NC}"
    fi

    # Run featureCounts for resistant samples of this pair
    if [ ${#RESIST_BAMS[@]} -gt 0 ]; then
        echo -e "${YELLOW}Counting resistant samples (${#RESIST_BAMS[@]} BAMs)${NC}"
        featureCounts -T $THREADS \
            -a "$GTF" \
            -o "$OUTDIR/rna_expression/${pair_id}_resist_counts.txt" \
            "${RESIST_BAMS[@]}"
    else
        echo -e "${RED}WARNING: No resistant BAMs found for pair $pair_id${NC}"
    fi

    # Create combined counts file for the pair (parent vs resist)
    if [ ${#PARENT_BAMS[@]} -gt 0 ] && [ ${#RESIST_BAMS[@]} -gt 0 ]; then
        echo -e "${YELLOW}Creating combined counts for pair $pair_id${NC}"
        featureCounts -T $THREADS \
            -a "$GTF" \
            -o "$OUTDIR/rna_expression/${pair_id}_combined_counts.txt" \
            -g gene_id \
            -t exon \
            --countReadPairs \
            -s 1 \
            "${PARENT_BAMS[@]}" "${RESIST_BAMS[@]}"
    fi
done

echo -e "${GREEN}Feature counting completed for all pairs${NC}"

# # # === DONE ===
# echo "\n[✓] Pipeline finished successfully. All outputs in: $OUTDIR"

# =======================
# DONE
# =======================
 end_time=$(date +%s)
# Calculate elapsed time
 elapsed=$((end_time - start_time))
# Convert to minutes and seconds
 minutes=$((elapsed / 60))
 seconds=$((elapsed % 60))
# echo "Pipeline completed in $ELAPSED seconds."
 echo "Pipeline completed in $minutes minute(s) and $seconds second(s)."