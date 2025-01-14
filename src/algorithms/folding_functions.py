import random as rd
from src.classes.protein import Protein
import src.visualisation.plot_functions as plot
from typing import Callable, Optional

def two_strings_fold(protein_sequence: str) -> list[int]:
    sequence_list: list[int] = []
    protein_length = len(protein_sequence)
    for i in range(protein_length):
        if i < (protein_length // 2 - 1):
            sequence_list.append(1)
        elif i == (protein_length // 2 - 1):
            sequence_list.append(2)
        elif i == (protein_length - 1):
            sequence_list.append(0)
        else:
            sequence_list.append(-1)
    return sequence_list

def random_fold(protein_sequence: str) -> list[int]:
    sequence_list: list[int] = []
    protein_length = len(protein_sequence)
    random_choice = [-2, -1, 1, 2]
    previous_direction = 0

    for _ in range(protein_length - 1):
        direction = rd.choice(random_choice)
        while direction == previous_direction:
            direction = rd.choice(random_choice)

        sequence_list.append(direction)
        previous_direction = -1 * direction

    sequence_list.append(0)
    return sequence_list

def random_iterated(sequence: str, fold_function: Callable[[str], list[int]], n: Optional[int] = 10000) -> None:
    score_list: list[int] = []
    best_structure: Optional[Protein] = None
    best_score: int = 0

    for _ in range(n):
        protein = Protein(sequence, fold_function)
        while protein.protein_rating == 1:
            protein = Protein(sequence, fold_function)
        score_list.append(protein.protein_rating)

        if protein.protein_rating < best_score:
            best_score = protein.protein_rating
            best_structure = protein
        
    plot.histogram(protein, score_list, save = True, iterations = n)

    if best_structure is not None:
        plot.visualize(best_structure, save = True)