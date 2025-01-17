from src.utils.helpers import save_and_visualize_results, random_fold
from src.classes.protein import Protein

class Random():
    """The Random class generates a random sequence for the folding direction"""

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Protein | None = None

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """Randomly generates sequences for a protein and calculates the scores."""
        for _ in range(repeats):
            self._random_iterated(iterations)

        # Save and visualize protein
        if self.best_protein is not None:
            save_and_visualize_results(self.best_protein, algorithm = "Random", histogram_data = self.histogram_data, 
            histogram = self.histogram_data[-1], iterations = iterations, show_plot= show_plot, save_plot= save_plot, save_data= save_data)
        else:
            print("Error: Did not find a valid protein.")

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