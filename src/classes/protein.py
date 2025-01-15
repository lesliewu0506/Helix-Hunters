import csv
from src.classes.grid import Grid
from src.classes.rating import Rating
from typing import Callable, Optional

class Protein():
    """
    A class to represent a protein and its attributes.
    It stores a protein sequence and uses a provided folding function
    to generate a folding pattern (amino_directions) or it uses a prebuilt folding pattern. Then it creates
    a structure (mapping of positions to amino acids) and computes
    the total rating of the protein.
    """

    def __init__(self, protein_sequence: str, amino_directions: Optional[list[int]] = None) -> None:
        self.protein_sequence: str = protein_sequence
        self.amino_directions: Optional[list[int]] = amino_directions
        self.structure: Optional[Grid] = None
        self.protein_rating: int = 1

    def build_structure(self, function: Callable[[str], list[int]]) -> None:
        """Creates the attributes for the protein with specific folding function."""
        self.amino_directions = function(self.protein_sequence)
        structure = Grid(self.protein_sequence, self.amino_directions)
        
        # Check if structure is valid else give rating 1
        if not structure.create_structure():
            self.protein_rating = 1
        else:
            self.structure = structure
            self.protein_rating = Rating(self.structure.get_structure()).get_rating()

    def build_no_function(self) -> None:
        """Builds the structure with given folding pattern."""
        if self.amino_directions is not None:
            structure = Grid(self.protein_sequence, self.amino_directions)

        # Check if structure is valid else give rating 1
        if not structure.create_structure():
            self.protein_rating = 1
        else:
            self.structure = structure
            self.protein_rating = Rating(self.structure.get_structure()).get_rating()

    def output_csv(self, file_path: Optional[str] = "output") -> None:
        """
        Creates a csv file containing the amino acids and their fold.
        Uses 'file_path' as output directory. 
        If not specified, creates 'output.csv' in current directory.
        """
        if self.amino_directions is not None:
            with open(f'{file_path}.csv', 'w', newline = '') as csvfile:
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