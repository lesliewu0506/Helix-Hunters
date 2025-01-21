from . import General
from src.utils import random_fold
from src.classes import Protein

class Random(General):
    """The Random class generates a random sequence for the folding direction"""

    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence)

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """Randomly generates sequences for a protein and calculates the scores."""
        self.run_algorithm(
            algorithm = "Random",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._random_iterated)

    def _random_iterated(self, n: int) -> None:
        """
        Iterates over multiple random generated folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(random_fold)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein.build_structure(random_fold)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score = protein.protein_rating
                self.best_protein = protein

        self.histogram_data.append(score_list)