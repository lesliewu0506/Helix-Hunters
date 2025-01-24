class Grid():
    """
    A class to represent a structure of a protein sequence.
    It contains a `dict` with:
    - keys (x, y, z): coordinates of an amino acid.
    - values (type, order): type of amino acid and the order in the chain.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).
    
    amino_directions : `list[int]`, optional
        A list of absolute directions.
    """
    # Maps direction to a change in x, y and z
    direction_map: dict[int, tuple[int, int, int]] = {
        0 : (0, 0, 0),
        1 : (1, 0, 0),
        -1 : (-1, 0, 0),
        2 : (0, 1, 0),
        -2 : (0, -1, 0),
        3 : (0, 0, 1),
        -3 : (0, 0, -1)}

    def __init__(self, protein_sequence: str, amino_directions: list[int] | None) -> None:
        """ 
        `structure` : `dict[tuple[int, int, int], tuple[str, int]]`
            The structure contains the coordinates of an amino acid and its order in the sequence.

        `protein_sequence` : str
            Protein sequence (for example `HHPHHHPH`).
    
        `amino_directions` : `list[int]`, optional
            A list of absolute directions.

        `x_current` : int
            The current x-coordinate. Default is `0`.

        `y_current` : int
            The current y-coordinate. Default is `0`.

        `z_current` : int
            The current z-coordinate. Default is `0`.
        """
        self.structure: dict[tuple[int, int, int], tuple[str, int]] = {}
        self.protein_sequence: str = protein_sequence
        self.amino_directions: list[int] | None = amino_directions

        self.x_current: int = 0
        self.y_current: int = 0
        self.z_current: int = 0

    def create_structure(self) -> bool:
        """
        Builds the dictionary structure.
        Returns True on success, else False.
        """
        if self.amino_directions is None:
            return False

        for i, amino in enumerate(self.protein_sequence):
            direction = self.amino_directions[i]
            # Check if amino added gives valid structure
            if not self._add_amino(amino, i, direction):
                return False
        return True
    
    def _add_amino(self, amino: str, i: int, direction: int) -> bool:
        """
        Adds amino acid to the structure.
        Return True on success, else False.
        """
        # Check if position is empty
        if (self.x_current, self.y_current, self.z_current) in self.structure:
            return False
        # Adds coordinates as keys with values [type of amino, order in chain]
        self.structure[(self.x_current, self.y_current, self.z_current)] = (amino, i)
        self._update_position(direction)
        return True

    def _update_position(self, direction: int) -> None:
        """Updates x, y and z based on direction."""
        dx, dy, dz = self.direction_map[direction]
        self.x_current += dx
        self.y_current += dy
        self.z_current += dz
    
    def get_structure(self) -> dict[tuple[int, int, int], tuple[str, int]]:
        """Returns the structure of the protein."""
        return self.structure