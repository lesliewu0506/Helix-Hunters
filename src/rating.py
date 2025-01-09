class Rating():
    """A class to represent the rating of a protein structure based on its sequence and structure."""

    def __init__(self, protein_sequence: str, structure: dict[tuple[int, int], tuple[str, int]]) -> None:
        self.protein_sequence: str = protein_sequence
        self.structure: dict[tuple[int, int], tuple[str, int]] = structure
        self.score: int = 0

        self._count_adjacent()

    def _count_adjacent(self) -> None:
        """Calculates the strength of the protein based on adjacent amino acids."""
        # Mapping for all neighbouring points
        direction_map = {(1, 0), (-1, 0), (0, 1), (0, -1)}

        for x_current, y_current in self.structure:

            amino_1 = self.structure[(x_current, y_current)][0]
            # Check for polar amino acid
            if amino_1 != 'P':
                # Find neighbouring coordinates
                for (dx, dy) in direction_map:
                    x_next = x_current + dx
                    y_next = y_current + dy

                    # Check for valid and non sequential points
                    if (x_next, y_next) in self.structure and not self._check_sequential(x_current, y_current, x_next, y_next):
                        amino_2 = self.structure[(x_next, y_next)][0]

                        self.score += self._check_pair(amino_1, amino_2)

        # Prevent double counting
        self.score = self.score // 2
    
    def _check_sequential(self, x_old: int, y_old: int, x_new: int, y_new: int) -> bool:
        """
        Checks if two amino acids are in a sequential order. 
        Returns True if sequential, False otherwise.
        """        
        i = self.structure[(x_old, y_old)][1]
        return self.structure[(x_new, y_new)][1] in [(i - 1), (i + 1)]

    def _check_pair(self, amino_1: str, amino_2: str) -> int:
        """Checks the pair for possible connections and returns the strength of connection."""
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return -1
        elif amino_1 == 'C' and amino_2 == 'C':
            return -5
        return 0

    def get_rating(self) -> int:
        """Returns rating of the protein."""
        return self.score