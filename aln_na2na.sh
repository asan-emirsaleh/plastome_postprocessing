#!/usr/bin/env bash

export PATH="${HOME}/miniconda3/envs/blast-env/bin:${PATH}"

# Argument parsing
#
# aa_seq=
# or
# aa_file = 
#
# na_ref =
# out_file = 
# out_format = (optional argument)

# aln
# translation table if needed = 11
# output to stdout if output file was not defined


usage() {
    echo "Usage: $0 [-s AA_SEQ | -f AA_FILE] -r NA_REF [-o OUT_FILE] [-t OUT_FORMAT]"
    exit 1
}

# Initialize variables
aa_seq=""
aa_file=""
na_ref=""
out_file=""
out_format=""

# Parse arguments
while getopts "s:f:r:o:t:" opt; do
    case $opt in
        s) aa_seq="$OPTARG" ;;
        f) aa_file="$OPTARG" ;;
        r) na_ref="$OPTARG" ;;
        o) out_file="$OPTARG" ;;
        t) out_format="$OPTARG" ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [ -z "$na_ref" ]; then
    echo "Error: Nucleotide reference (-r) is required"
    usage
fi

if [ -z "$aa_seq" ] && [ -z "$aa_file" ]; then
    echo "Error: Either amino acid sequence (-s) or file (-f) must be provided"
    usage
fi

if [ -n "$aa_seq" ] && [ -n "$aa_file" ]; then
    echo "Error: Cannot specify both sequence and file"
    usage
fi

# Prepare input
if [ -n "$aa_seq" ]; then
    echo "$aa_seq" > temp_aa.fasta
    input_file="temp_aa.fasta"
else
    input_file="$aa_file"
fi

# Prepare output redirection
output_redirect=""
if [ -n "$out_file" ]; then
    output_redirect="-o $out_file"
fi

# Create BLAST database if needed
if [ ! -f "${na_ref}.nhr" ]; then
    makeblastdb -in "$na_ref" -dbtype nucl
fi

# Run alignment
blastn \
    -query "$input_file" \
    -db "$na_ref" \
    -outfmt 0 \
    -word_size 7 \
    ${out_file:+-out "$out_file"}

# Cleanup if temporary file was created
if [ -n "$aa_seq" ]; then
    rm temp_aa.fasta
fi