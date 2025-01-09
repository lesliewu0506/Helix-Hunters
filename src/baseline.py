from protein import Protein
import plot_functions as plot
import folding_functions as fold

def main(sequence, fold_function):
    protein = Protein(sequence, fold_function)
    print(protein.protein_rating)
    plot.visualize(protein)
    protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    main(protein_sequence, fold.two_strings_fold)