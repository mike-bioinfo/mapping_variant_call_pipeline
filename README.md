# # README
Use chmod +xr on the bash script, make sure to create "results" directory and modify all the path required by the script to run successfully.
The script also can run without rna-seq replicates directly to WES or WGS data

Current Script Behavior
The script will:

✅ Process WES/WGS data normally (alignment, deduplication, BQSR)

✅ Perform variant calling (Mutect2) on WES/WGS data

✅ Skip RNA-seq alignment when rna_reps is empty/NA

✅ Skip featureCounts for samples with no RNA BAM files

Usage:

$./mapping_DNA-RNA_variant-call.sh
