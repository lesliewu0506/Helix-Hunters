import itertools
import multiprocessing
import csv

import src.visualisation.plot_functions as plot
from src.classes.protein import Protein

def brute_force(sequence: str, save: bool = False) -> None:
    """
    Function that brute forces every possible combination.
    Plots the best structure and prints the rating of the best structure.
    Finally saves the directions into a csv file.
    """
    # num_processes: int = multiprocessing.cpu_count()
    num_processes: int = 8
    best_score = 1
    # Load data
    foldings: list[list[int]] = read_csv(sequence)

    # Use multiprocessing pool to distribute work load
    with multiprocessing.Pool(num_processes) as pool:
        results = pool.map(evaluate_folding_wrapper, [(sequence, folding) for folding in foldings])

    # Process results
    best_structures: list[Protein] = []
    i = 0
    for rating, structure in results:
        if rating < best_score:
            best_score = rating
            protein = structure
            best_structures = []
        elif rating == best_score and protein not in best_structures:
            best_structures.append(protein)
             
    # Plot best structure and save data
    for i, protein in enumerate(best_structures):
        plot.visualize(protein, show = False, save = save, file_path = f"src/brute_force/protein_structures/{sequence}_{i}")
        protein.output_csv(file_path = f"src/brute_force/best_folding/{sequence}_{i}")
    print(f"Best rating: {best_score}")

def read_csv(protein_sequence: str) -> list[list[int]]:
    """Reads csv file and returns the directions as list in list."""
    foldings: list[list[int]] = []

    with open(f"src/brute_force/all_foldings/{protein_sequence}.csv", 'r') as csvfile:
        reader = csv.reader(csvfile)

        for row in reader:
            foldings.append([int(x) for x in row])

    return foldings
    
def evaluate_folding_wrapper(args: tuple[str, list[int]]) -> tuple[int, Protein]:
    """Wrapper function to unpack arguments for evaluate_folding."""
    sequence, folding = args
    return evaluate_folding(sequence, folding)

def evaluate_folding(sequence: str, folding: list[int]) -> tuple[int, Protein]:
    """Evaluate a single folding and return the result."""
    protein = Protein(sequence, amino_directions = folding)
    rating = protein.get_rating()
    return rating, protein