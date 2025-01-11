import itertools
import multiprocessing
import csv
import pandas as pd

import plot_functions as plot
from protein import Protein
from grid import Grid
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

def refine_csv(sequence: str) -> None:
    """Filter the csv by removing rows with invalid prefixes."""
    df = pd.read_csv(f"{sequence}.csv", header = None)
    # Split df into equal sizes
    cpu_cores = multiprocessing.cpu_count()
    chunk_size = len(df) // cpu_cores
    chunks = [df.iloc[i : i + chunk_size] for i in range(0, len(df), chunk_size)]

    # Use multiprocessing pool to distribute work load
    with multiprocessing.Pool(cpu_cores) as pool:
        results = pool.map(_process_chunk, chunks)
    
    # Retrieve results
    valid_rows = (itertools.chain.from_iterable(result for result in results))
    pd.DataFrame(valid_rows).to_csv(f"{sequence}_refined.csv", index = False, header = False)

def _process_chunk(chunk: pd.DataFrame) -> list[list[int]]:
    """Process a chunk of the dataframe and returns valid rows."""
    invalid_prefixes: set[tuple[int]] = set()
    valid_rows: list[list[int]] = []

    # Iterate over all rows with index
    for _, row in chunk.iterrows():
        # Convert to list
        row_list = row.tolist()
        
        # Check if sequence contains invalid prefix
        if any(tuple(row_list[:len(prefix)]) == prefix for prefix in invalid_prefixes):
            continue

        # Check if sequence has new invalid prefix
        invalid_prefix = _check_valid_folding(row_list)
        if invalid_prefix is not None:
            invalid_prefixes.add(tuple(invalid_prefix))
        else:
            valid_rows.append(row_list)

    return valid_rows

def _check_valid_folding(folding: list[int]) -> Optional[list[int]]:
    """
    Helper function that checks a folding sequence.
    If folding is not valid, return the invalid prefix.
    Else return None.
    """
    dummy_sequence = 'H' * len(folding)
    grid = Grid(dummy_sequence, folding)

    if not grid.create_structure():
        return grid.invalid_prefix

    return None
    
def evaluate_folding_wrapper(args: tuple[str, list[int]]) -> tuple[int, Protein]:
    """Wrapper function to unpack arguments for evaluate_folding."""
    sequence, folding = args
    return evaluate_folding(sequence, folding)

def evaluate_folding(sequence: str, folding: list[int]) -> tuple[int, Protein]:
    """Evaluate a single folding and return the result."""
    protein = Protein(sequence, amino_directions = folding)
    rating = protein.get_rating()
    return rating, protein