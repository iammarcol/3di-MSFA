# this script is used to modify the output 3di MSA, since the DF_Colores is made
# based on the positions in the reference sequence, so this script basically
# adjust 3di MSA and this adjusted MSA is the input for making DF_Colores_(3di)/MODIF table 

import sys

def process_multifasta(file_path):
    sequences = {}
    reference_sequence = None
    current_sequence_id = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_sequence_id = line[1:]
                sequences[current_sequence_id] = ''
                if reference_sequence is None:
                    reference_sequence = current_sequence_id
            else:
                sequences[current_sequence_id] += line

    return sequences, reference_sequence

def process_reference_sequence(reference_sequence):
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    positions = [i for i, char in enumerate(reference_sequence) if char.upper() in amino_acids]
    return positions

def modify_sequences(sequences, reference_positions):
    modified_sequences = {}

    for sequence_id, sequence in sequences.items():
        modified_sequence = ''.join(sequence[pos] if pos < len(sequence) else '-' for pos in reference_positions)
        modified_sequences[sequence_id] = modified_sequence

    return modified_sequences

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py input_file.fasta")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = "MSA_3di_reference.fasta"

    sequences, reference_sequence_id = process_multifasta(input_file)

    if reference_sequence_id not in sequences:
        print(f"Reference sequence {reference_sequence_id} not found.")
        return

    reference_positions = process_reference_sequence(sequences[reference_sequence_id])
    modified_sequences = modify_sequences(sequences, reference_positions)

    with open(output_file, 'w') as out_file:
        for sequence_id, modified_sequence in modified_sequences.items():
            out_file.write(f">{sequence_id}\n{modified_sequence}\n")

    print(f"MSA modified sequences saved to {output_file}")

if __name__ == "__main__":
    main()
