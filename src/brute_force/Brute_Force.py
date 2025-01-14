import itertools
import multiprocessing
import csv

import src.visualisation.plot_functions as plot
from src.classes.protein import Protein
from src.classes.grid import Grid
from typing import Optional

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

def generate_all_foldings(protein_sequence: str) -> None:
    """
    Generate all possible foldings for a given sequence length 
    where the first item is always 1 and the last item 0
    and consecutive items are never opposing directions.
    Then saves it to a csv file.
    """
    sequence_length = len(protein_sequence)
    directions = [0, 1, 2]

    with open(f'{protein_sequence}.csv', 'w', newline = '') as csvfile:
        
        writer = csv.writer(csvfile)

        for folding in itertools.product(directions, repeat = sequence_length - 2):
            result = _check_folding(folding)
            if result is not None:
                writer.writerow(result)

    csvfile.close()

def _check_folding(folding: list[int]) -> Optional[list[int]]:
    """
    Checks if a configuration is valid.
    Translates the folding list first into absolute directions.
    Returns the list of directions if valid, else None.
    """
    abs_folding: list[int] = _direction_translator(folding)
    invalid_prefix: Optional[list[int]] = _check_valid_folding(abs_folding)

    if invalid_prefix is None:
        return abs_folding
    return None

def _direction_translator(directions: list[int]) -> list[int]:
    """
    Helper function for translating the relative paths (0, 1, 2),
    to absolute paths (-2, -1, 1, 2). Returns list of directions.
    """
    direction_map: dict[int, list[int]] = {1: [2, 1, -2], -1: [-2, -1, 2], 2: [-1, 2, 1], -2: [1, -2, -1]}
    folding_sequence: list[int] = [1]

    for i, direction in enumerate(directions):
        folding_sequence.append(direction_map[folding_sequence[i]][direction])

    folding_sequence.append(0)
    return folding_sequence

def _check_valid_folding(folding: list[int]) -> Optional[list[int]]:
    """
    Helper function that checks a folding sequence.
    If folding is not valid, returns the invalid prefix.
    Else return None.
    """
    dummy_sequence = 'H' * len(folding)
    grid = Grid(dummy_sequence, folding)

    if grid.create_structure():
        return None

    return grid.invalid_prefix
    
def evaluate_folding_wrapper(args: tuple[str, list[int]]) -> tuple[int, Protein]:
    """Wrapper function to unpack arguments for evaluate_folding."""
    sequence, folding = args
    return evaluate_folding(sequence, folding)

def evaluate_folding(sequence: str, folding: list[int]) -> tuple[int, Protein]:
    """Evaluate a single folding and return the result."""
    protein = Protein(sequence, amino_directions = folding)
    rating = protein.get_rating()
    return rating, protein