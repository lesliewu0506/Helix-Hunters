class Rating():
    """
    A class to represent the rating of a protein structure
    based on its sequence and structure.
    It will take in a protein structure
    and calculate its rating based on rules:
    - If two `H` amino acids are next to each other, it gets a score `-1`.
    - If `H` and `C` are next to each other, it also gets a score `-1`.
    - If two `C` amino acids are next to each other, its gets a score `-5`.
    
    This only applies if the amino acids are not sequential.
    
    Parameter
    ---------
    structure : `dict[tuple[int, int, int], tuple[str, int]]`
        The structure contains the information of the folding of a protein.
    """

    def __init__(self, structure: dict[tuple[int, int, int], tuple[str, int]]) -> None:
        """
        Attributes
        ----------
        `structure` : `dict[tuple[int, int, int], tuple[str, int]]`
            The structure contains the coordinates
            of an amino acid and its order in the sequence.

        `score` : int
            The score based on its structure.
        """
        self.structure: dict[tuple[int, int, int], tuple[str, int]] = structure
        self.score: int = 0

        self._count_adjacent()

    def _count_adjacent(self) -> None:
        """
        Calculates the strength of the protein based on adjacent amino acids.
        """
        # Mapping for all neighbouring points
        direction_map = {(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)}

        for (x_current, y_current, z_current), (amino_1, order_1) in self.structure.items():

            # Skip if polar amino acid
            if amino_1 == 'P':
                continue

            # Find neighbouring coordinates
            for (dx, dy, dz) in direction_map:
                x_next = x_current + dx
                y_next = y_current + dy
                z_next = z_current + dz

                # Search pairs
                pair = self.structure.get((x_next, y_next, z_next))
                if pair is not None:
                    amino_2, order_2 = pair

                    # Skip sequential points
                    if order_2 not in (order_1 - 1, order_1 + 1):
                        self.score += self._check_pair(amino_1, amino_2)

        # Prevent double counting
        self.score = self.score // 2

    def _check_pair(self, amino_1: str, amino_2: str) -> int:
        """
        Checks the pair for possible connections
        and returns the strength of connection.
        """
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return -1
        elif amino_1 == 'C' and amino_2 == 'C':
            return -5
        return 0

    def get_rating(self) -> int:
        """Returns rating of the protein."""
        return self.score