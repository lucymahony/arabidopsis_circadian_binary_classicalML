#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 1 
#SBATCH -c 16 
#SBATCH --mem 6G			
#SBATCH --time=0-10:00				
#SBATCH --output output_%x		
#SBATCH --mail-type=END,FAIL			
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	
#SBATCH -p ei-medium

source package b0ed0698-358b-4c9b-9d21-603ea8d6e478  # bedtools 2.31.0
source package 9117dc0d-5bb0-48ce-83c3-bc0db08e0cd3  # gffread 0.12.8 

# Input files
DIR="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data"
genes_labels_file="$DIR/raw/pnas.2103070118.sd01.csv"
genome_assembly_file="$DIR/raw/GCF_000001735.3_TAIR10_genomic.fna"
genome_annotation_file="$DIR/raw/genomic.gff"
araport_annotation_file="$DIR/raw/Araport11_GFF3_genes_transposons.current.gff"

# Output files
transcript_ids="$DIR/processed/transcripts.txt"
filtered_annotation="$DIR/processed/filtered.gtf"
transcript_fasta="$DIR/processed/transcripts.fa"
promoter_bed="$DIR/processed/promoters.bed"
promoter_fasta="$DIR/processed/promoters.fa"
found_ids="$DIR/processed/found_transcripts.txt"
missing_ids="$DIR/processed/missing_transcripts.txt"


#  Extract transcript IDs (uppercase, remove duplicates)
tail -n +2 "$genes_labels_file" | cut -d',' -f1 | awk '{ print toupper($1) }' | sort -u > "$transcript_ids"
feature_list="antisense_lncRNA antisense_RNA lnc_RNA miRNA miRNA_primary_transcript mRNA ncRNA pseudogenic_transcript pseudogenic_tRNA rRNA snoRNA snRNA tRNA"
LC_ALL=C awk -v idfile="$transcript_ids" '
BEGIN {
    FS = "\t"; OFS = "\t";
    # Build the set of allowed feature types
    keep["antisense_lncRNA"]
    keep["antisense_RNA"]
    keep["lnc_RNA"]
    keep["miRNA"]
    keep["miRNA_primary_transcript"]
    keep["mRNA"]
    keep["ncRNA"]
    keep["pseudogenic_transcript"]
    keep["pseudogenic_tRNA"]
    keep["rRNA"]
    keep["snoRNA"]
    keep["snRNA"]
    keep["tRNA"]
    while ((getline < idfile) > 0) ids[$1] = 1;
}
($3 in keep) {
    match($9, /ID=([^;]+)/, a);
    tid = a[1];
    if (tid in ids) print;
}' "$araport_annotation_file" > "$filtered_annotation"


# check for missing IDs in the filtered annotation 
LC_ALL=C awk -v idfile="$transcript_ids" '
BEGIN {
    FS = "\t"; OFS = "\t";
    while ((getline < idfile) > 0) ids[$1] = 1;
}
{
    match($9, /ID=([^;]+)/, a);
    tid = a[1];
    if (tid in ids) {
        delete ids[tid]
    }
}
END {
    for (id in ids) {
        print id
    }
}' "$filtered_annotation" > "$missing_ids"

# The missing IDs are difficult ones often mitochondrial genes - remove transcript number and search both id and in annotation file. 
# For the missing IDs check between ID= and . 
difficult_filtered="$DIR/processed/difficult_filtered.gtf"
# remove any .1 or .2 from missing IDs 
sed -i 's/\.[0-9]*//g' "$missing_ids"
LC_ALL=C awk -v idfile="$missing_ids" '
BEGIN {
    FS = "\t"; OFS = "\t";
    keep["antisense_lncRNA"]
    keep["antisense_RNA"]
    keep["lnc_RNA"]
    keep["miRNA"]
    keep["miRNA_primary_transcript"]
    keep["mRNA"]
    keep["ncRNA"]
    keep["pseudogenic_transcript"]
    keep["pseudogenic_tRNA"]
    keep["rRNA"]
    keep["snoRNA"]
    keep["snRNA"]
    keep["tRNA"]
    keep["uORF"] 
    while ((getline < idfile) > 0) ids[$1] = 1;
}
($3 in keep) {
    match($9, /ID=([^.]+)/, a)
    tid = a[1];
    if (tid in ids) print;
}' "$araport_annotation_file" > "$difficult_filtered"

# Append the difficult filtered to the filtered annotation file
cat "$difficult_filtered" >> "$filtered_annotation"
# Print the number of rows in the filtered annotation file
echo "Number of rows in filtered annotation file: $(wc -l < "$filtered_annotation")"
# Convert the genome assembly file chromosome names to Chr1 e.c.t NC_003070.9
filtered_annotation_renamed="$DIR/processed/filtered_renamed.gtf"
LC_ALL=C awk 'BEGIN {
    map["Chr1"] = "NC_003070.9"
    map["Chr2"] = "NC_003071.7"
    map["Chr3"] = "NC_003074.8"
    map["Chr4"] = "NC_003075.7"
    map["Chr5"] = "NC_003076.8"
    map["ChrC"] = "NC_001284.2"
    map["ChrM"] = "NC_000932.1"
    FS = OFS = "\t"
}
{
    if ($1 in map) $1 = map[$1]
    print
}' $filtered_annotation > $filtered_annotation_renamed

gffread "$filtered_annotation_renamed" -g "$genome_assembly_file" -w "$transcript_fasta"
LC_ALL=C awk -v idfile="$transcript_ids" '
BEGIN {
    FS = "\t"; OFS = "\t";
    feature_list = "antisense_lncRNA antisense_RNA lnc_RNA miRNA miRNA_primary_transcript mRNA ncRNA pseudogenic_transcript pseudogenic_tRNA rRNA snoRNA snRNA tRNA"
    n = split(feature_list, features, " ")
    for (i = 1; i <= n; i++) keep[features[i]] = 1
    while ((getline < idfile) > 0) {
        raw = toupper($1)
        ids[raw] = 1
        fallback[raw ".1"] = 1
    }
}
($3 in keep) {
    match($9, /ID=([^;]+)/, a)
    tid = toupper(a[1])
    if (!(tid in ids || tid in fallback)) next
    chr = $1
    strand = $7
    if (strand == "+") {
        start = $4 - 1500
        end = $4 - 1 
    } else {
        start = $5
        end = $5 + 1499
    }
    if (start < 1) start = 1
    print chr, start - 1, end, tid, ".", strand
}' "$araport_annotation_file" > "$promoter_bed"

promoter_bed_renamed="$DIR/processed/promoter_bed_renamed.gtf"
LC_ALL=C awk 'BEGIN {
    map["Chr1"] = "NC_003070.9"
    map["Chr2"] = "NC_003071.7"
    map["Chr3"] = "NC_003074.8"
    map["Chr4"] = "NC_003075.7"
    map["Chr5"] = "NC_003076.8"
    map["ChrC"] = "NC_001284.2"
    map["ChrM"] = "NC_000932.1"
    FS = OFS = "\t"
}
{
    if ($1 in map) $1 = map[$1]
    print
}' $promoter_bed > $promoter_bed_renamed

# Extract promoter sequences
bedtools getfasta -fi "$genome_assembly_file" -bed "$promoter_bed_renamed" -s -name -fo "$promoter_fasta"

echo "Done. Outputs:"
echo " - Transcript FASTA: $transcript_fasta"
echo " - Promoter FASTA:   $promoter_fasta"