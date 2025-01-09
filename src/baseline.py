import plot_functions as plot
import folding_functions as fold
from protein import Protein
from typing import Callable

def main(sequence: str, fold_function: Callable[[str], list[int]]) -> None:
    protein = Protein(sequence, fold_function)
    
    while protein.protein_rating == 1:
        protein = Protein(sequence, fold_function)

    plot.visualize(protein)
    # protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    main(protein_sequence, fold.random_fold)