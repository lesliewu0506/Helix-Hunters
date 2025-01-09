from protein import Protein
import plot_functions

def two_strings_fold(protein_sequence):
    sequence_list = []
    protein_length = len(protein_sequence)
    for i in range(protein_length):
        if i < (protein_length // 2 - 1):
            sequence_list.append(1)
        elif i == (protein_length // 2 - 1):
            sequence_list.append(2)
        elif i == (protein_length - 1):
            sequence_list.append(0)
        else:
            sequence_list.append(-1)
    return sequence_list

def main(sequence, fold_function):
    protein = Protein(sequence, fold_function)
    print(protein.protein_rating)
    # plot = Plot(protein)
    # plot.plot_structure()
    protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    main(protein_sequence, two_strings_fold)