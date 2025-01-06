def rating():
    return

def two_strings_fold(protein_sequence):
    sequence_list = []
    protein_length = len(protein_sequence)
    for i in range(protein_length):
        if i < (protein_length / 2):
            sequence_list.append(1)
        elif i == (protein_length / 2):
            sequence_list.append(2)
        else:
            sequence_list.append(-1)
    return sequence_list

if __name__ == "__main__":
    Protein = "HHPHHHPH"
    sequence = two_strings_fold(Protein)
    for i in range(len(Protein)):
        print(f"{Protein[i], sequence[i]}")
