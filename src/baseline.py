import plot_functions as plot
import folding_functions as fold

from protein import Protein
from Brute_Force import brute_force, generate_all_foldings, _check_valid_folding, refine_csv
from typing import Callable, Optional

def random_iterated(sequence: str, fold_function: Callable[[str], list[int]]) -> None:
    score_list: list[int] = []
    best_structure: Optional[Protein] = None
    best_score: int = 0

    for _ in range(1000):
        protein = Protein(sequence, fold_function)
        while protein.protein_rating == 1:
            protein = Protein(sequence, fold_function)
        score_list.append(protein.protein_rating)

        if protein.protein_rating < best_score:
            best_score = protein.protein_rating
            best_structure = protein
        
    plot.histogram(score_list)

    if best_structure is not None:
        plot.visualize(best_structure)

def main(sequence: str, fold_function: Callable[[str], list[int]]) -> None:
    protein = Protein(sequence, fold_function)

    while protein.protein_rating == 1:
        protein = Protein(sequence, fold_function)
    plot.visualize(protein)
    # protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HHPHHHPHPH"
    generate_all_foldings(protein_sequence)
    refine_csv(protein_sequence)
    # brute_force(protein_sequence)
    # protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    # main(protein_sequence, fold.random_fold)
    # random_iterated(protein_sequence, fold.random_fold)