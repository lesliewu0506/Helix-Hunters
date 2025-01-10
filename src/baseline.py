import itertools

import plot_functions as plot
import folding_functions as fold
from protein import Protein

def generate_all_foldings(sequence_length: int) -> list[list[int]]:
    """Generate all possible foldings for a given sequence length where the last item is always 0."""
    directions = [-2, -1, 1, 2]
    all_foldings = [[1] + list(folding) + [0] for folding in itertools.product(directions, repeat=sequence_length - 2)]
    print(all_foldings[0])
    return all_foldings

def main() -> None:
    pass

if __name__ == "__main__":
    protein_sequence = "HHPHHHPHPHHHPH"
    print(len(generate_all_foldings(14)))