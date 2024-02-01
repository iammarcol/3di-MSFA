#!/bin/bash

# check for arguments num
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 PDB_DIR MSA_FILE"
    exit 1
fi

# pdb files directory
PDB_DIR="$1"

# 3d descriptors output dir
OUTPUT_DIR="$PDB_DIR/descriptors"

DI_DIR="$PDB_DIR/3di_seq"

# create output dir 
mkdir -p "$OUTPUT_DIR"
mkdir -p "$PDB_DIR/3di_seq"

# iterate over PDB files
for pdb_file in "$PDB_DIR"/*.pdb; do
    # extract PDB filename without extension
    pdb_filename=$(basename "$pdb_file")

    # run FoldSeek on each PDB file and save the output in the output directory
    # for this to work foldseek has to be configured
    foldseek structureto3didescriptor "$pdb_file" "$OUTPUT_DIR/${pdb_filename}_3di.dump"
done

for dump_file in "$OUTPUT_DIR"/*.dump; do
    dump_filename=$(basename "$dump_file" .dump)
    awk -F'\t' '{print ">"$1"\n"$3; next}{print}' "$dump_file" > "$DI_DIR/${dump_filename}_seq.fasta"
done

# here I make multiple 3di fastafile

for fasta_file in "$DI_DIR"/*.fasta; do
    # create a temporary file with modified fasta ids
    temp_file="$DI_DIR/temp.fasta"
    sed 's/.pdb$//' "$fasta_file" > "$temp_file"

    # concatenate the modified fasta file to the output file
    cat "$temp_file" >> "$DI_DIR/3di_multiplefasta.fasta"

    # cemove the temporary file
    rm "$temp_file"
done

msa_3di="$DI_DIR/3di_multiplefasta.fasta"

echo "Finished generating multiple 3di fasta file! This fasta doesn't have gaps."
echo "Making MSA 3di based on the original MSA..."

# if you want the alignment to be generated as well unlock these lines
# clustalo -i "$DI_DIR/3di_multiplefasta.fasta" -o "$DI_DIR/3di_msa.fasta"
# echo "Finished generating 3di alignment file!"

# Python script
msa_file="$2"
output_file="./MSA_3di.fasta"
python aa_msa_to_3di.py "$msa_file" "$msa_3di" "$output_file"
python mod_mfasta.py "$output_file"