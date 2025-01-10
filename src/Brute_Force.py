import itertools
import multiprocessing

import plot_functions as plot
from protein import Protein

def brute_force(sequence: str) -> None:
    num_processes: int = multiprocessing.cpu_count()
    foldings = generate_all_foldings(len(sequence))
    best_score = 1

    # Use multiprocessing pool and map to distribute work load
    with multiprocessing.Pool(num_processes) as pool:
        results = pool.map(evaluate_folding_wrapper, [(sequence, folding) for folding in foldings])

    # Process results
    for rating, structure in results:
        if rating < best_score:
            best_score = rating
            protein = structure

    plot.visualize(protein)
    print(f"Best rating: {best_score}")
    protein.output_csv(file_name = f"best_fold_{sequence}")

def generate_all_foldings(sequence_length: int) -> list[list[int]]:
    """Generate all possible foldings for a given sequence length where the last item is always 0,
    and consecutive items are never opposing directions."""
    directions = [-2, -1, 1, 2]  # Possible folding directions excluding 0
    all_foldings = []

    for folding in itertools.product(directions, repeat = sequence_length - 2):
        valid = True
        for i in range(1, len(folding)):
            if folding[i] == -folding[i - 1]:  # Check if consecutive items are opposing
                valid = False
                break
        if valid:
            all_foldings.append([1] + list(folding) + [0])
    print(len(all_foldings))
    return all_foldings

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