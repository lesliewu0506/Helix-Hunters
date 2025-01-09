import plot_functions as plot
import folding_functions as fold
from protein import Protein
from typing import Callable

def main(sequence: str, fold_function: Callable[[str], list[int]]) -> None:
    protein = Protein(sequence, fold_function)
    print(protein.protein_rating)
    plot.visualize(protein)
    protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    main(protein_sequence, fold.two_strings_fold)