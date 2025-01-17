import itertools
import csv

from src.utils.helpers import direction_translator
from src.classes.grid import Grid
from typing import Optional

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
    invalid_prefix: Optional[list[int]] = _check_valid_folding(abs_folding)

    if invalid_prefix is None:
        return abs_folding
    return None

def _check_valid_folding(folding: list[int]) -> Optional[list[int]]:
    """
    Helper function that checks a folding sequence.
    If folding is not valid, returns [0].
    Else return None.
    """
    dummy_sequence = 'H' * len(folding)
    grid = Grid(dummy_sequence, folding)

    if grid.create_structure():
        return None

    return [0]