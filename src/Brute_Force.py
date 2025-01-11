import itertools
import multiprocessing
import csv

import plot_functions as plot
from protein import Protein
from typing import Optional

def brute_force(sequence: str, foldings: list[list[int]], save: Optional[bool] = False) -> None:
    """
    Function that brute forces every possible combination.
    Plots the best structure and prints the rating of the best structure.
    Finally saves the directions into a csv file.
    """
    num_processes: int = multiprocessing.cpu_count()
    best_score = 1

    # Use multiprocessing pool to distribute work load
    with multiprocessing.Pool(num_processes) as pool:
        results = pool.map(evaluate_folding_wrapper, [(sequence, folding) for folding in foldings])

    # Process results
    for rating, structure in results:
        if rating < best_score:
            best_score = rating
            protein = structure

    # Plot best structure and save data
    plot.visualize(protein, save)
    print(f"Best rating: {best_score}")
    protein.output_csv(file_name = f"best_fold_{sequence}")

def generate_all_foldings(protein_sequence: str) -> None:
    """
    Generate all possible foldings for a given sequence length 
    where the first item is always 1 and the last item 0
    and consecutive items are never opposing directions.
    Then saves it to a csv file.
    """
    sequence_length = len(protein_sequence)
    directions = [-2, -1, 1, 2]

    with open(f'{protein_sequence}.csv', 'w', newline = '') as csvfile:
        
        writer = csv.writer(csvfile)

        for folding in itertools.product(directions, repeat = sequence_length - 2):
            result = _check_valid_sequence(folding)
            if result is not None:
                writer.writerow(result)

    csvfile.close()

def _check_valid_sequence(folding: list[int]) -> Optional[list[int]]:
    """
    Helper function for checking valid sequence.
    Returns the complete directions if valid, else None.
    """
    # Checks if the first one is valid
    if folding[0] == -1:
        return None
    
    # Check if consecutive items are opposing
    for i in range(1, len(folding)):
        if folding[i] == -folding[i - 1]:
            return None

    return ([1] + list(folding) + [0])

def evaluate_folding_wrapper(args: tuple[str, list[int]]) -> tuple[int, Protein]:
    """Wrapper function to unpack arguments for evaluate_folding."""
    sequence, folding = args
    return evaluate_folding(sequence, folding)

def evaluate_folding(sequence: str, folding: list[int]) -> tuple[int, Protein]:
    """Evaluate a single folding and return the result."""
    protein = Protein(sequence, amino_directions = folding)
    rating = protein.get_rating()
    return rating, protein