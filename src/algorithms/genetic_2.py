import random as rd

from src.utils.helpers import random_fold
from src.classes.protein import Protein
from typing import Callable, Optional
from .general import General

class Genetic_Algorithm():
    """The Genetic random class generates a sequence for the folding direction."""

    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence)

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """The Genetic algorithm generates sequences for a protein and calculates the scores"""
        self.run_algorithm(
            algorithm = "Genetic",
            show_plot = show_plot,
            save_plot = save_plot,
            saved_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._greedy_iterated)
        
    def ted(self, n: int) -> None:
s



initialized_population()
evaluate_fit()
select_parents()
generate_crossover()
apply_mutation()
check_convergence()