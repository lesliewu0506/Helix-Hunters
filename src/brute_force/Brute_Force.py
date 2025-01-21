import multiprocessing
import csv

from src.visualisation import plot
from src.classes import Protein

protein_sequence_map: dict[str, str] = {
    "HHPHHHPHPHHHPH" : "1",
    "HPHPPHHPHPPHPHHPPHPH" : "2",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

def brute_force(sequence: str, save: bool = False) -> None:
    """
    Function that brute forces every possible combination.
    Plots the best structure and prints the rating of the best structure.
    Finally saves the directions into a csv file.
    """
    folder = protein_sequence_map[sequence]
    num_processes: int = multiprocessing.cpu_count()
    best_score = 0
    best_structures: list[Protein] = []

    # Open CSV once
    with open(f"src/brute_force/all_foldings/{sequence}.csv", 'r') as csvfile:
        reader = csv.reader(csvfile)
        chunk: list[list[int]] = []

        for row in reader:
            directions = [int(x) for x in row]
            chunk.append(directions)

            # Process chunk when chunk size reached
            if len(chunk) == 1000000:
                best_score, best_structures = process_chunk(sequence, chunk, best_score, best_structures, num_processes)
                chunk = []

        # Process final chunk
        if chunk:
            best_score, best_structures = process_chunk(sequence, chunk, best_score, best_structures, num_processes)

    csvfile.close()

    # Plot best structure and save data
    for i, protein in enumerate(best_structures):
        plot.visualize_protein(protein, "Brute Force", show = False, save = save, file_path = f"src/brute_force/protein_structures/{folder}/{sequence}_{i}")
        protein.output_csv(file_path = f"src/brute_force/best_folding/{folder}/{sequence}_{i}")

    print(f"Best rating: {best_score}")

def process_chunk(sequence: str, foldings: list[list[int]], best_score: int, best_structures: list[Protein], num_processes: int) -> tuple[int, list[Protein]]:
    """Processes chunk in parallel and updates best_score and best_structures."""
    with multiprocessing.Pool(num_processes) as pool:
        results = pool.map(evaluate_folding_wrapper, [(sequence, f) for f in foldings])

    for rating, protein in results:
        if rating < best_score:
            best_score = rating
            best_structures = [protein]

        elif rating == best_score:
            best_structures.append(protein)

    return best_score, best_structures

def evaluate_folding_wrapper(args: tuple[str, list[int]]) -> tuple[int, Protein]:
    """Wrapper function to unpack arguments for evaluate_folding."""
    sequence, folding = args
    return evaluate_folding(sequence, folding)

def evaluate_folding(sequence: str, folding: list[int]) -> tuple[int, Protein]:
    """Evaluate a single folding and return the result."""
    protein = Protein(sequence, amino_directions = folding)
    protein.build_no_function()
    rating = protein.get_rating()
    return rating, protein