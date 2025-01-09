class Grid():
    """
    A class to represent a structure of a protein sequence.
    It contains self.structure, a dict with:
    - keys (x, y): coordinates of an amino acid.
    - values (type, order): type of amino acid and the order in the chain.
    """

    def __init__(self, protein_sequence: str, amino_directions: list[int]) -> None:
        self.structure: dict[tuple[int, int], tuple[str, int]] = {}
        self.protein_sequence: str = protein_sequence
        self.amino_directions: list[int] = amino_directions
        self.x_current: int = 0 
        self.y_current: int = 0

    def create_structure(self) -> bool:
        """
        Builds the dictionary structure.
        Returns True on success, else False.
        """
        for i in range(len(self.protein_sequence)):
            if not self._add_amino(self.protein_sequence[i], i, self.amino_directions[i]):
                return False
        return True
    
    def _add_amino(self, amino: str, i: int, direction: int) -> bool:
        """
        Adds amino acid to the structure.
        Return True on success, else False.
        """
        # Check if position is empty
        if not self._check_fold(self.x_current, self.y_current):
            return False
        # Adds coordinates as keys with values [type of amino, order in chain]
        self.structure[(self.x_current, self.y_current)] = (amino, i)
        self._update_position(direction)
        return True

    def _update_position(self, direction: int) -> None:
        """Updates x and y based on direction."""
        # Maps direction to a change in x and y
        direction_map = {0 : (0, 0), 1 : (1, 0), -1 : (-1, 0), 2: (0, 1), -2: (0, -1)}

        dx, dy = direction_map[direction]
        self.x_current += dx
        self.y_current += dy

    def _check_fold(self, x_new: int, y_new: int) -> bool:
        """
        Checks if the folding structure is valid.
        Return True on success, else False.
        """
        return not (x_new, y_new) in self.structure

    def get_structure(self) -> dict[tuple[int, int], tuple[str, int]]:
        """Returns the structure of the protein."""
        return self.structure