import csv
from grid import Grid
from rating import Rating
from typing import Callable

class Protein():
    """
    A class to represent a protein and its attributes.
    It stores a protein sequence and uses a provided folding function
    to generate a folding pattern (amino_directions). Then it creates
    a structure (mapping of positions to amino acids) and computes
    the total rating of the protein.
    """

    def __init__(self, protein_sequence: str, function: Callable[[str], list[int]]) -> None:
        self.protein_sequence: str = protein_sequence
        self.amino_directions: list[int] = []
        self.structure: dict[tuple[int, int], tuple[str, int]] = {}
        self.protein_rating: int = 0
        self._build_structure(function)

    def _build_structure(self, function: Callable[[str], list[int]]) -> None:
        """Creates the attributes for the protein with specific fold."""
        self.amino_directions = function(self.protein_sequence)

        self.structure = Grid(self.protein_sequence, self.amino_directions).get_structure()
        self.protein_rating = Rating(self.protein_sequence, self.structure).get_rating()

    def output_csv(self) -> None:
        """Creates a csv file containing the amino acids and their fold."""
        with open('output.csv', 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            # Mandatory header line
            writer.writerow(['amino', 'fold'])

            for amino, direction in zip(self.protein_sequence, self.amino_directions):
                writer.writerow([amino, direction])

            # Mandatory footer line
            writer.writerow(['score', self.protein_rating])
    
    def get_rating(self) -> int:
        """Return the rating of the protein fold."""
        return self.protein_rating