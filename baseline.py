import src.visualisation.plot_functions as plot
import src.algorithm.folding_functions as fold

from src.classes.protein import Protein
from src.brute_force.Brute_Force import brute_force, generate_all_foldings
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
        plot.visualize(best_structure, save = True)

def main(sequence: str) -> None:
    generate_all_foldings(sequence)
    print("Done generating")
    brute_force(sequence, save = True)

if __name__ == "__main__":
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    # main(protein_sequence)
    random_iterated(protein_sequence, fold.random_fold)