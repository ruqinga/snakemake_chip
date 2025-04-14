#!/bin/bash

# Usage: nohup bash download_sra.sh <metadata_file> &
# The <metadata_file> should have the following format:
#
# study_accession    run_accession    sample_title
# PRJNA592545        SRR10560107      E13.5_female_H3K27me3_rep1
# PRJNA592545        SRR10560113      E13.5_female_H3K27me3_rep2
# PRJNA592545        SRR13296473      E13.5_female_H3K27me3_rep3


# Check if metadata argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <metadata_file>"
    exit 1
fi

metadata="$1"  # The metadata file is passed as the first argument
metadata_dir=$(dirname "$metadata_file")

# Change to the directory containing the metadata file
cd "$metadata_dir" || { echo "Failed to change directory to $metadata_dir"; exit 1; }

# Extract SRA list
awk 'NR > 1 {print $2}' "${metadata}" > sra_list.txt

# Download SRA files
iseq -i sra_list.txt -g > iseq.log 2>&1

# Check if files are downloaded
missing_files=false
while IFS=$'\t' read -r study_accession run_accession sample_title; do
    # Check for expected files
    for file in "${run_accession}"*.fastq.gz; do
        if [ ! -f "$file" ]; then
            echo "File $file not found!"
            missing_files=true
        fi
    done
done < "${metadata}"

# List files before renaming
ls -lh ./ >> sra_bytes.log

# Rename files
while IFS=$'\t' read -r study_accession run_accession sample_title; do
    for file in "${run_accession}"*.fastq.gz; do
        if [ -f "$file" ]; then
            # Rename file
            new_name=$(echo "$file" | sed "s/^${run_accession}/${sample_title}/")
            mv "$file" "$new_name"
            echo "Renamed $file to $new_name"
        else
            echo "File $file not found, skipping..."
        fi
    done
done < "${metadata}"

# List files after renaming
ls -lh ./ >> renamed_bytes.log