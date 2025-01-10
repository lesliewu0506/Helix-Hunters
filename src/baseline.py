import itertools
import multiprocessing

import plot_functions as plot
import folding_functions as fold
from protein import Protein

def generate_all_foldings(sequence_length: int) -> list[list[int]]:
    """Generate all possible foldings for a given sequence length where the last item is always 0."""
    directions = [-2, -1, 1, 2]
    all_foldings = [[1] + list(folding) + [0] for folding in itertools.product(directions, repeat=sequence_length - 2)]
    print(len(all_foldings))
    return all_foldings

def evaluate_folding(sequence: str, folding: list[int], result_queue: multiprocessing.Queue) -> None:
    """Evaluate a single folding and put the result in the queue."""
    protein = Protein(sequence, folding)
    rating = protein.get_rating()
    return rating, protein.structure

def main(sequence: str, num_processes: int) -> None:
    pass

if __name__ == "__main__":
    protein_sequence = "HHPH"
    main(protein_sequence, num_processes = multiprocessing.cpu_count())