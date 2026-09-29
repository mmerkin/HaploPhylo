#!/bin/bash
# Haplophylo v1.3
# Usage: bash Haplophylo.sh -r reference.fa -a annotation.bed -i /path/to/bams -o output_prefix -f depth.tsv -n tree_build -s 20 
set -euo pipefail

# Set variables

log_file="logfile_haplophylo.txt"
> "$log_file"
BQ=20
MQ=20
N_THRESHOLD=20
reference=""
annotation=""
input_path=""
output_prefix=""
sample_file=""
min_depth=""
min_allele_depth=""
allele_balance=""


# Send EVERYTHING (stdout + stderr) to log file
exec >> "$log_file" 2>&1

# Functions

show_help() {
  echo "Usage: $0 -r <reference> -a <annotation> -i <input_path> -o <output_prefix> -f <sample_file> [-n snp_call,consensus,tree_build] [-s 10] [-d 3] [-b 0.2]"
  echo
  echo "Options:"
  echo "  -r    Mitochondrial reference genome (required)"
  echo "  -a    Annotation bed file (required)"
  echo "  -i    Path to bam files (required)"
  echo "  -o    Output prefix (required)"
  echo "  -f    File containing sample names (required) and depths (optional)"
  echo "  -s    Minimum depth to keep sample (optional)"
  echo "  -d    Minimum mean depth for genotype call (optional)"
  echo "  -b    Maximum allele balance allowed for genotype call (optional)"
  echo "  -n    Skip a time-consuming step in the pipeline"
  echo "  -h    Show help message"
}

log() {
    local msg="$1"
    local dest="${2:-both}"

    case "$dest" in
        terminal)
            echo -e "$msg" > /dev/tty
            ;;
        log)
            echo -e "$msg" >> "$log_file"
            ;;
        both)
            echo -e "$msg" > /dev/tty
            echo -e "$msg" >> "$log_file"
            ;;
        *)
            echo "Invalid log destination: $dest" >&2
            ;;
    esac
}

print_banner() {
clear > /dev/tty
cat << 'EOF' > /dev/tty
          _   _    _    ____  _     ___  ____  _   ___   ___     ___  
         | | | |  / \  |  _ \| |   / _ \|  _ \| | | \ \ / / |   / _ \ 
         | |_| | / _ \ | |_) | |  | | | | |_) | |_| |\ V /| |  | | | |
         |  _  |/ ___ \|  __/| |__| |_| |  __/|  _  | | | | |__| |_| |
         |_| |_/_/   \_\_|   |_____\___/|_|   |_| |_| |_| |_____\___/ 
                                                                  
EOF
}


# Arguments

declare -A skip

while getopts "r:a:i:o:f:s:d:b:n:h" opt; do
  case $opt in
    r) reference="$OPTARG" ;;
    a) annotation="$OPTARG" ;;
    i) input_path="$OPTARG" ;;
    o) output_prefix="$OPTARG" ;;
    f) sample_file="$OPTARG" ;;
    s) min_depth="$OPTARG" ;;
    d) min_allele_depth="$OPTARG" ;;
    b) allele_balance="$OPTARG" ;;
    n)
      IFS=',' read -ra steps <<< "$OPTARG"
      for step in "${steps[@]}"; do
        case "$step" in
          snp_call|consensus|tree_build|variant) # Options that can be skipped
            skip["$step"]=1
            ;;
          *)
            log "Unknown step to skip: $step"
            show_help > /dev/tty
            exit 1
            ;;
        esac
      done
      ;;

    h)
      show_help > /dev/tty
      exit 0
      ;;
    \?)
      log "Invalid option: -$OPTARG"
      show_help > /dev/tty
      exit 1
      ;;
    :)
      log "Option -$OPTARG requires an argument"
      show_help > /dev/tty
      exit 1
      ;;
  esac
done

trap 'log "Command failed: $BASH_COMMAND" terminal' ERR


if [[ -z "$reference" || -z "$annotation" || -z "$input_path" || -z "$output_prefix" || -z "$sample_file" ]]; then
  log "ERROR: Missing required arguments"
  show_help > /dev/tty
  exit 1
fi

missing=()
for cmd in bcftools bedtools samtools iqtree snp-sites; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        missing+=("$cmd")
    fi
done

if ((${#missing[@]} > 0)); then
    log "There are missing commands. Have you activated the conda environment?"
    exit 1
fi


run_cmd() {
    local cmd="$*"
    echo "$cmd" >> "$log_file"
    eval "$cmd"
}

# Suppress annoying warning

bcftools() {
    command bcftools "$@" \
        2> >(grep -vF '[W::bcf_hdr_check_sanity] MQ should be declared as Type=Float' >&2)
}


# Start

print_banner

log "Script started at $(date)" log
log "Starting run of Haplophylo\nPlease report any errors as an issue on github\n"
log "Reference: $reference"
log "Annotation file: $annotation"
log "Input path: $input_path"
log "Output prefix: $output_prefix"
log "Sample file: $sample_file"
[[ -n "$min_depth" ]] && log "Minimum sample depth: $min_depth"
[[ -n "$min_allele_depth" ]] && log "Minimum allele depth: $min_allele_depth"
[[ -n "$allele_balance" ]] && log "Maximum allele balance: $allele_balance"
log "\nBegin analysis"

# Find samples

mkdir -p ${output_prefix}_consensus
mkdir -p ${output_prefix}_variants

cols=$(awk 'NR>1 {print NF; exit}' "$sample_file")

if [[ -n "$min_depth" && "$cols" -lt 2 ]]; then
  log "ERROR: The sample file only contains one column. A 2-column tsv file is required to filter by sample depth with the -s option"
  exit 1
fi

if [[ -z "$min_depth" ]]; then
  mapfile -t SAMPLES < <(awk 'NR>1 {print $1}' "$sample_file")
else
  mapfile -t SAMPLES < <(awk -v min_dp="$min_depth" 'NR>1 && $2 > min_dp {print $1}' "$sample_file")
fi

retained=${#SAMPLES[@]}
total=$(awk 'NR>1' "$sample_file" | wc -l)

log "Using $retained samples out of $total"

for sample in "${SAMPLES[@]}"; do
    sample_base=$(basename "$sample" .bam)
    printf "sample_base='%s'\n" "$sample_base"
done

# Call variants


BAMS=""
for sample in "${SAMPLES[@]}"; do
    bam="$input_path/$sample"
    if [ -f "$bam" ]; then
        BAMS="$BAMS $bam"
    fi
done

if [ -z "$BAMS" ]; then
    log "ERROR: No bam files available for variant calling"
    exit 1
fi

for BAM in $BAMS; do
sample=$(basename "$BAM" .bam)
echo -e "$BAM\t$sample"
done > ${output_prefix}_variants/${output_prefix}_samples.txt

if [[ ! -v 'skip[snp_call]' ]]; then
log "Calling SNPs"

run_cmd bcftools mpileup \
-f "$reference" \
-q $MQ -Q $BQ \
--ignore-RG \
-a AD,DP \
-Ou \
$BAMS | \
bcftools call \
--ploidy 1 \
-mv \
-Oz \
-o "${output_prefix}_variants/${output_prefix}_variants.raw.vcf.gz"

run_cmd bcftools index ${output_prefix}_variants/${output_prefix}_variants.raw.vcf.gz

# Rename samples

run_cmd bcftools query -l ${output_prefix}_variants/${output_prefix}_variants.raw.vcf.gz |
while read bam; do
basename "$bam" .bam
done > ${output_prefix}_variants/${output_prefix}_new_names.txt

run_cmd bcftools reheader \
-s ${output_prefix}_variants/${output_prefix}_new_names.txt \
-o ${output_prefix}_variants/${output_prefix}_variants.renamed.vcf.gz \
${output_prefix}_variants/${output_prefix}_variants.raw.vcf.gz

# Remove indels
run_cmd bcftools view -v snps -Oz -o ${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz ${output_prefix}_variants/${output_prefix}_variants.renamed.vcf.gz
run_cmd bcftools index ${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz

snp_count=$(bcftools view -H ${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz | wc -l)
log "Called $snp_count SNPs"

else
log "Skipping SNP calling"
fi

# Make gene counts output file

run_cmd bedtools intersect -a "$annotation" -b ${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz -c > ${output_prefix}_gene_variant_counts.tsv

# Create depth per variant summary file for masking

create_depthpersite_tsv() {
local vcf="$1"
local sample_base="$2"
local depthpersite_tsv="$3"
{
echo -e "CHROM\tPOS\tDepth(A)\tDepth(C)\tDepth(G)\tDepth(T)\tcalled_allele\tallele_balance"
bcftools query \
-s "$sample_base" \
-f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%GT]\n' \
"$vcf" |
awk -F'\t' -v OFS='\t' '{
chrom = $1
pos = $2
ref = $3
alt = $4
ad = $5
gt = $6

# Initialise depths
A = 0
C = 0
G = 0
T = 0

# Number of ALT alleles
nalt = split(alt, alts, ",")

# REF is first in AD
nalleles = nalt + 1

# Build allele list
allele[1] = ref

for (i = 1; i <= nalt; i++)
allele[i + 1] = alts[i]

# Split AD
nad = split(ad, depths, ",")

# Assign depths to A/C/G/T
for (i = 1; i <= nalleles; i++) {
if (i > nad)
continue
depth = depths[i]
if (depth == ".")
depth = 0

if (allele[i] == "A")
A = depth
else if (allele[i] == "C")
C = depth
else if (allele[i] == "G")
G = depth
else if (allele[i] == "T")
T = depth
}

# Determine called allele from haploid GT
if (gt == "0")
called_allele = ref
else if (gt ~ /^[1-9][0-9]*$/) {
gt_index = gt + 1
called_allele = allele[gt_index]
}
else
called_allele = "."

# Calculate total depth
total_depth = A + C + G + T

# Calculate allele balance
if (total_depth > 0) {
if (called_allele == "A")
called_depth = A
else if (called_allele == "C")
called_depth = C
else if (called_allele == "G")
called_depth = G
else if (called_allele == "T")
called_depth = T
else
called_depth = 0

allele_balance = called_depth / total_depth
}
else
allele_balance = "1"

print chrom, pos, A, C, G, T, called_allele, allele_balance
}'
} > "$depthpersite_tsv"
}


total_masked=0
processed_samples=0
# Create summary file for masking
printf 'Sample\tnMiss\tnLC\tnAB\tTotal\n' > "${output_prefix}_filtering_summary.tsv"

for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
depthpersite_tsv="${output_prefix}_consensus/${sample_base}_depthpersite.tsv"
snp_vcf="${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz"

create_depthpersite_tsv \
"$snp_vcf" \
"$sample_base" \
"$depthpersite_tsv"

# Define masks

MASK="${output_prefix}_consensus/${sample_base}_mask.bed"
NO_DATA_MASK="${output_prefix}_consensus/${sample_base}_no_data.bed"
LOW_DP_MASK="${output_prefix}_consensus/${sample_base}_low_depth.bed"
AB_MASK="${output_prefix}_consensus/${sample_base}_allele_balance.bed"
> "$MASK"
> "$NO_DATA_MASK"
> "$LOW_DP_MASK"
> "$AB_MASK"

# Mask positions that are missing

awk 'NR > 1 && $3 == 0 && $4 == 0 && $5 == 0 && $6 == 0 {
print $1 "\t" ($2 - 1) "\t" $2
}' "$depthpersite_tsv" >> "$NO_DATA_MASK"

n_miss=$(wc -l < "$NO_DATA_MASK")

# Mask positions with low depth

if [[ -n "${min_allele_depth:-}" ]]; then
awk -v min_dp="$min_allele_depth" '
NR > 1 {
total_dp = $3 + $4 + $5 + $6
if (total_dp > 0 && total_dp < min_dp)
print $1 "\t" ($2-1) "\t" $2
}' "$depthpersite_tsv" >> "$LOW_DP_MASK"
fi

n_lowC=$(wc -l < "$LOW_DP_MASK")

# Mask positions with low allele balance

if [[ -n "${allele_balance:-}" ]]; then
awk -v min_ab="$allele_balance" '
NR > 1 {
if ($8 < min_ab)
print $1 "\t" ($2-1) "\t" $2
}' "$depthpersite_tsv" >> "$AB_MASK"
fi

n_lowAB=$(wc -l < "$AB_MASK")

# Merge masks

sort -k1,1 -k2,2n -u \
"$NO_DATA_MASK" "$LOW_DP_MASK" "$AB_MASK" > "$MASK"
masked_sites=$(wc -l < "$MASK")
echo "Sample $sample_base masked sites: $masked_sites"
total_masked=$((total_masked + masked_sites))
processed_samples=$((processed_samples + 1))
printf '%s\t%d\t%d\t%d\t%d\n' "$sample_base" "$n_miss" "$n_lowC" "$n_lowAB" "$masked_sites" >> "${output_prefix}_filtering_summary.tsv"

done

if [[ "$processed_samples" -gt 0 ]]; then
avg_masked=$(awk -v total="$total_masked" -v n="$processed_samples" 'BEGIN {printf "%.2f", total/n}')
log "Average masked sites per sample: $avg_masked"
else
log "No samples processed"
fi

# Create consensus sequence and extract gene sequences

if [[ ! -v 'skip[consensus]' ]]; then
log "Generating consensus sequences"

for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
bcftools consensus \
-f "$reference" \
-s "$sample_base" \
-m "${output_prefix}_consensus/${sample_base}_mask.bed" \
${output_prefix}_variants/${output_prefix}_variants.snps.vcf.gz \
> "${output_prefix}_consensus/${sample_base}_full_mt_consensus.fasta" 2> >(grep -v 'Note: applying IUPAC codes' >&2)


rm -f "${output_prefix}_consensus/${sample_base}_full_mt_consensus.fasta.fai"
bedtools getfasta -fi "${output_prefix}_consensus/${sample_base}_full_mt_consensus.fasta" \
-bed "$annotation" -nameOnly -fo - | \
awk -v sample="$sample_base" 'BEGIN{RS=">"; ORS=""} NR>1 {
split($0, lines, "\n")
gene = lines[1]
seq = ""
for(i=2;i<=length(lines);i++) seq = seq lines[i]
print ">" sample "|" gene "\n" seq "\n"
}' > "${output_prefix}_consensus/${sample_base}_genes_consensus.fasta"
done

else
log "Skipping generation of consensus sequences"
fi

# Create parition file and concatenated gene file


concat_sequence="${output_prefix}_all_sequences.fa"
partition_file="${output_prefix}_partitions.nex"
> "$concat_sequence"
> "$partition_file"
pos=1

first_sample="${SAMPLES[0]}"  # Use first sample to get gene lengths
first_sample=$(basename "$first_sample" .bam)

# Read gene names from BED

gene_list=()
while read -r chrom start end gene; do
gene_list+=("$gene")
done < "$annotation"

# Initialize per-sample concatenated strings

declare -A concat_seqs
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
concat_seqs[$sample_base]=""
done

# Loop through genes
for gene in "${gene_list[@]}"; do
seq=$(awk -v sample="$first_sample" -v gene="$gene" '
BEGIN {found=0; seq=""}
/^>/ {
gsub(/^ +| +$/,"",$0)
if ($0 ~ "^>" sample "\\|" gene "$") {found=1} else {found=0}
next
}
found {seq=seq $0}
END {print seq}' "${output_prefix}_consensus/${first_sample}_genes_consensus.fasta")
seq_len=$(echo -n "$seq" | tr -d '\n' | wc -c)

if [ "$seq_len" -eq 0 ]; then
log "Warning: gene $gene sequence length is 0 in sample $first_sample"
exit 1
fi

end_pos=$((pos + seq_len - 1))
echo "DNA, $gene = $pos-$end_pos" >> "$partition_file"

# Append sequences for all samples
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
gene_seq=$(awk -v sample_base="$sample_base" -v gene="$gene" '
BEGIN {found=0; seq=""}
/^>/ {
gsub(/^ +| +$/,"",$0)
if ($0 ~ "^>" sample_base "\\|" gene "$") {found=1} else {found=0}
next
}
found {seq=seq $0}
END {print seq}' "${output_prefix}_consensus/${sample_base}_genes_consensus.fasta")
concat_seqs[$sample_base]="${concat_seqs[$sample_base]}$gene_seq"
done
pos=$((end_pos + 1))
done

# Write concatenated fasta
for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
echo ">$sample_base" >> "$concat_sequence"
echo "${concat_seqs[$sample_base]}" >> "$concat_sequence"
done

# Final sanity check

ALN_LEN=$(awk '/^[^>]/ {print length($0); exit}' "$concat_sequence")
LAST_PART=$(tail -n1 "$partition_file" | awk -F"=" '{print $2}' | tr -d ' ' | awk -F"-" '{print $2}')

if [ "$ALN_LEN" -ne "$LAST_PART" ]; then
    log "ERROR: alignment length ($ALN_LEN) does not match last partition end ($LAST_PART). Please report this bug."
    exit 1
fi

GENOME_LEN=$(awk '!/^>/ {total += length($0)} END {print total}' "$reference")
PERCENTAGE=$(awk -v u="$ALN_LEN" -v g="$GENOME_LEN" 'BEGIN {printf "%.2f", (u/g)*100}')

log "\nAlignment length: $ALN_LEN (${PERCENTAGE}% of genome used)"

# snp sites only

run_cmd snp-sites "$concat_sequence" > ${output_prefix}_variants_only_alignment.fa

# Make new vcf for concat_sequence

MAP="${output_prefix}_variants/${output_prefix}_coordinate_map.tsv"
GENEVCF="${output_prefix}_variants/${output_prefix}_gene_only_snps.vcf"

awk 'BEGIN {
OFS="\t"
print "GENE", "CHROM", "BED_START", "BED_END", \
"SEQ_LENGTH", "CONCAT_START", "CONCAT_END"
offset = 0
}
{
chrom = $1
bed_start = $2
bed_end = $3
gene = $4
seq_length = bed_end - bed_start
concat_start = offset + 1
concat_end = offset + seq_length
print gene, chrom, bed_start, bed_end, \
seq_length, concat_start, concat_end
offset += seq_length
}' "$annotation" > "$MAP"

TOTAL_LENGTH=$(awk '
NR > 1 { total += $5 }
END { print total }
' "$MAP")

zcat $snp_vcf | awk \
-v map="$MAP" \
-v total_length="$TOTAL_LENGTH" '

BEGIN {
OFS="\t"
snp_count = 0

# Read coordinate map
while ((getline line < map) > 0) {
if (line ~ /^GENE\t/)
continue

split(line, a, "\t")
gene = a[1]
gchrom[gene] = a[2]
gstart[gene] = a[3]
gend[gene] = a[4]

cstart[gene] = a[6]
n++
genes[n] = gene
}
close(map)
}
/^#/ {
if ($0 ~ /^#CHROM/) {
print "##contig=<ID=concatenated,length=" total_length ">"
print "##INFO=<ID=ORIG_CHROM,Number=1,Type=String,Description=\"Original VCF chromosome\">"
print "##INFO=<ID=ORIG_POS,Number=1,Type=Integer,Description=\"Original VCF position\">"
print "##INFO=<ID=GENE,Number=1,Type=String,Description=\"Gene corresponding to concatenated coordinate\">"
}
print
next
}

# VCF records

{
chrom = $1
pos = $2

# Check every gene because genes can overlap
for (i = 1; i <= n; i++) {
gene = genes[i]
# BED [start,end) corresponds to VCF positions
# start+1 through end
if (chrom == gchrom[gene] &&
pos > gstart[gene] &&
pos <= gend[gene]) {
# Convert genomic coordinate to concatenated coordinate
newpos = cstart[gene] + \
(pos - gstart[gene] - 1)
# Add original coordinate information
info = $8
if (info == ".")
info = ""
if (info != "")
info = info ";"
info = info \
"ORIG_CHROM=" chrom \
";ORIG_POS=" pos \
";GENE=" gene
# Increment SNP counter
snp_count++
$1 = "concatenated"
$2 = newpos
$3 = "SNP" snp_count
$8 = info
print
}
}
}
' > "$GENEVCF"
run_cmd bgzip -f "$GENEVCF"
run_cmd bcftools index "${GENEVCF}.gz"

for sample in "${SAMPLES[@]}"; do
sample_base=$(basename "$sample" .bam)
depthpersite_tsv="${output_prefix}_consensus/${sample_base}_gene_only_depthpersite.tsv"
create_depthpersite_tsv \
"${GENEVCF}.gz" \
"$sample_base" \
"$depthpersite_tsv"
done

# Check Ns

OUTFILE="${output_prefix}_high_missing.tsv"

: > "$OUTFILE"

read AVG_N COUNT <<EOF
$(awk -v thresh="$N_THRESHOLD" -v out="$OUTFILE" '
BEGIN {
  OFS="\t"
  print "sample_id", "N_count" > out
}

/^>/ {
  if (seqs > 0) {
    total += n
    if (n > thresh) {
      print id, n >> out
      count++
    }
  }
  id = substr($0,2)
  seqs++
  n = 0
  next
}

{
  n += gsub(/[Nn]/, "&")
}

END {
  total += n
  if (n > thresh) {
    print id, n >> out
    count++
  }

  if (seqs > 0)
    print total / seqs, count
  else
    print 0, 0
}
' "$concat_sequence")
EOF

log "Average number of Ns: $AVG_N"

if (( COUNT > 0 )); then
  log "WARNING: $COUNT samples found with more than $N_THRESHOLD missing sites."
  log "See $OUTFILE for details"
fi


# IQTree

if [[ ! -v 'skip[tree_build]' ]]; then
echo "Creating phylogenetic tree"

mkdir -p ${output_prefix}_tree

run_cmd iqtree3 -s "$concat_sequence" -p "$partition_file" -m MFP -bb 1000 -T AUTO -pre ${output_prefix}_tree/${output_prefix}_tree > /dev/null 2>&1

echo "Check ${output_prefix}_tree/${output_prefix}_tree.log for iqtree progress" 

cp ${output_prefix}_tree/${output_prefix}_tree.treefile ${output_prefix}_tree.treefile

log "Run finished successfully. Find output tree in ${output_prefix}_tree.treefile"

else
log "Skipping creation of phylogenetic tree"
fi

log "Script finished at $(date)" log

# End of script
