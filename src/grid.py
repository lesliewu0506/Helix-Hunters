class Grid():
    """
    A class to represent a structure of a protein sequence.
    It contains self.structure, a dict with:
    - keys (x, y): coordinates of an amino acid.
    - values (type, order): type of amino acid and the order in the chain.
    """

    def __init__(self, protein_sequence: str, amino_directions: list[int]) -> None:
        self.structure: dict[tuple[int, int], tuple[str, int]] = self._create_structure(protein_sequence, amino_directions)

    def _create_structure(self, protein_sequence: str, amino_directions: list[int]) -> dict[tuple[int, int], tuple[str, int]]:
        """Builds the dictionary structure and returns the structure."""
        structure: dict[tuple[int, int], tuple[str, int]] = {}    
        x_current= 0 
        y_current= 0

        for i in range(len(protein_sequence)):
            # Adds coordinates as keys with values [type of amino, order in chain]
            structure[(x_current, y_current)] = (protein_sequence[i], i)
            x_current, y_current = self._update_position(x_current, y_current, amino_directions[i])

        return structure
    
    def _update_position(self, x_old: int, y_old: int, direction: int) -> tuple[int, int]:
        """Updates x and y based on direction."""
        # Maps direction to a change in x and y
        direction_map = {0 : (0, 0), 1 : (1, 0), -1 : (-1, 0), 2: (0, 1), -2: (0, -1)}

        dx, dy = direction_map[direction]
        return x_old + dx, y_old + dy

    def get_structure(self) -> dict[tuple[int, int], tuple[str, int]]:
        """Returns the structure of the protein."""
        return self.structure