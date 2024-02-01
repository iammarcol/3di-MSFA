#!/usr/bin/env python3
import sys

def read_fasta(file_path):
    sequences = {}
    current_id = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ''
            else:
                sequences[current_id] += line
    return sequences

def match_sequences_order(msa_sequences, di3_sequences):
    matched_sequences = {}
    for seq_id_msa, amino_acid_sequence in msa_sequences.items():
        for seq_id_di3, di3_sequence in di3_sequences.items():
            if seq_id_msa.lower() in seq_id_di3.lower() or seq_id_di3.lower() in seq_id_msa.lower():
                matched_sequences[seq_id_msa] = di3_sequence
                break
    return matched_sequences

def replace_amino_acids_with_3di(msa_sequences, di3_sequences, ids_to_process):
    modified_sequences = {}
    for seq_id in ids_to_process:
        amino_acid_sequence_msa = msa_sequences.get(seq_id, '')
        di3_sequence = di3_sequences.get(seq_id, '')
        modified_sequence = ''

        i_msa, i_di3 = 0, 0

        while i_msa < len(amino_acid_sequence_msa) and i_di3 < len(di3_sequence):
            aa_msa = amino_acid_sequence_msa[i_msa]
            aa_di3 = di3_sequence[i_di3]

            if aa_msa == '-':
                modified_sequence += '-'
                i_msa += 1
            elif aa_di3 == '-':
                modified_sequence += '-'
                i_di3 += 1
            else:
                modified_sequence += aa_di3
                i_msa += 1
                i_di3 += 1

        modified_sequences[seq_id] = modified_sequence

    return modified_sequences

def adjust_sequence_length(sequences):
    max_length = max(len(seq) for seq in sequences.values())
    adjusted_sequences = {seq_id: seq.ljust(max_length, '-') for seq_id, seq in sequences.items()}
    return adjusted_sequences

def write_fasta(output_file, msa_sequences, adjusted_sequences):
    with open(output_file, 'w') as file:
        for seq_id, sequence in msa_sequences.items():
            modified_sequence = adjusted_sequences.get(seq_id, '')
            file.write(f'>{seq_id}\n{modified_sequence}\n')

if len(sys.argv) != 4:
    print("Usage: python script.py input_file1 input_file2 output")
    sys.exit(1)

msa_file_path = sys.argv[1]
di3_file_path = sys.argv[2]
output_file_path = sys.argv[3]

# Read sequences from files
msa_sequences = read_fasta(msa_file_path)
di3_sequences = read_fasta(di3_file_path)

# Match sequences order based on IDs
matched_sequences = match_sequences_order(msa_sequences, di3_sequences)

# Get IDs from di3_sequences to process
ids_to_process = set(di3_sequences.keys())

# Replace amino acids with 3di characters, preserving gaps, for selected IDs
modified_sequences = replace_amino_acids_with_3di(msa_sequences, matched_sequences, ids_to_process)

# Adjust sequences to the same length
adjusted_sequences = adjust_sequence_length(modified_sequences)

# Write the modified sequences to a new file while maintaining the order of msa_sequences
write_fasta(output_file_path, msa_sequences, adjusted_sequences)

print(f'Modified MSA written to {output_file_path}')
