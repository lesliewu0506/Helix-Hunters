import itertools
import csv

from src.utils import direction_translator
from typing import Optional

# Maps direction to a change in x and y
direction_map = {
    0 : (0, 0),
    1 : (1, 0),
    -1 : (-1, 0),
    2: (0, 1),
    -2: (0, -1)
}

def generate_all_foldings(protein_sequence: str) -> None:
    """
    Generate all possible foldings for a given sequence length 
    where the first item is always 1 and the last item 0
    and consecutive items are never opposing directions.
    Then saves it to a csv file.
    """
    sequence_length = len(protein_sequence)
    directions = [0, 1, 2]

    with open(f'{protein_sequence}.csv', 'w', newline = '') as csvfile:
        
        writer = csv.writer(csvfile)

        for folding in itertools.product(directions, repeat = sequence_length - 2):
            result = _check_folding(list(folding))
            if result is not None:
                writer.writerow(result)

    csvfile.close()

def _check_folding(folding: list[int]) -> Optional[list[int]]:
    """
    Checks if a configuration is valid.
    Translates the folding list first into absolute directions.
    Returns the list of directions if valid, else None.
    """
    abs_folding: list[int] = direction_translator(folding)

    if _check_valid_folding(abs_folding):
        return abs_folding
    return None

def _check_valid_folding(folding: list[int]) -> bool:
    """
    Helper function that checks a folding sequence.
    If folding is not valid, returns False.
    Else return True.
    """
    coordinates: set[tuple[int, int]] = set()
    x_current: int = 0
    y_current: int = 0

    # Add coordinates
    for direction in folding:
        if (x_current, y_current) in coordinates:
            return False
        coordinates.add((x_current, y_current))
        
        # Update positions
        dx, dy = direction_map[direction]
        x_current += dx
        y_current += dy
    return True