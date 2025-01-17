import random as rd

from src.utils.helpers import save_and_visualize_results
from src.classes.protein import Protein

class Greedy():
    """The Greedy random class generates a sequence for the folding direction."""

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Protein | None = None

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """Greedily generates sequences for a protein and calculates the scores."""
        for _ in range(repeats):
            self._greedy_iterated(iterations)

        # Save and visualize protein
        save_and_visualize_results(self.best_protein, algorithm = "Greedy", histogram_data = self.histogram_data, 
        histogram = self.histogram_data[-1], iterations = iterations, show_plot= show_plot, save_plot= save_plot, save_data= save_data)

    def _greedy_iterated(self, n: int) -> None:
        """
        Iterates over multiple random generated folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(self._greedy_fold)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein = Protein(self.protein_sequence)
                protein.build_structure(self._greedy_fold)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score = protein.protein_rating
                self.best_protein = protein

        self.histogram_data.append(score_list)

    def _greedy_fold(self, protein_sequence: str) -> list[int]:
        """
        Generates a random folding sequence.
        Uses relative directions (0, 1, 2) and translates them 
        into absolute directions (-2, -1, 1, 2).
        Returns the sequence as a list.
        """
        relative_direction_list: list[int] = []
        directions = [0, 1, 2]

        for i in range(len(protein_sequence) - 2):
            
            if i % 5 == 0:
                direction = self._check_greedy(relative_direction_list)
                relative_direction_list.append(direction)

            else:
                direction = rd.choice(directions)
                relative_direction_list.append(direction)

        return self._direction_translator(relative_direction_list)
    
    def _direction_translator(self, directions: list[int]) -> list[int]:
        """
        Helper function for translating the relative paths (0, 1, 2),
        to absolute paths (-2, -1, 1, 2). Returns list of directions.
        """
        direction_map: dict[int, list[int]] = {1: [2, 1, -2], -1: [-2, -1, 2], 2: [-1, 2, 1], -2: [1, -2, -1]}
        folding_sequence: list[int] = [1]

        for i, direction in enumerate(directions):
            folding_sequence.append(direction_map[folding_sequence[i]][direction])

        folding_sequence.append(0)
        return folding_sequence

    def _check_greedy(self, direction_list: list[int]) -> int:
        """
        Implements a greedy algorithm to find the best possible next direction.
        Evaluates all possible continuations and selecting lowest score directions randomly if there are ties.
        """
        # Evaluate all possible next directions (0, 1, 2) and collect scores
        protein_scores = []
        best_directions = []
        minimum_score = 1

        for direction in [0, 1, 2]:
            new_direction_list = direction_list + [direction]
            abs_direction = self._direction_translator(new_direction_list)
            abs_direction.remove(0)

            # Create a Protein object and compute its rating
            protein = Protein(self.protein_sequence[:len(new_direction_list)], amino_directions=abs_direction)
            score = protein.protein_rating
            protein_scores.append(score)

            # Update best directions based on the current score
            if score < minimum_score:
                minimum_score = score
                best_directions = [direction]

            elif score == minimum_score:
                best_directions.append(direction)

        # Randomly choose one of the best directions
        return rd.choice(best_directions)