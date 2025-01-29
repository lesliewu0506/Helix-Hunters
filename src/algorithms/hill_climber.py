import random as rd

from . import General
from src.classes import Protein
from src.utils import (random_fold, MAX_UNCHANGED_ITERATIONS,
                       DIRECTION_CHOICES_2D, DIRECTION_CHOICES_3D)
from typing import Callable

class HillClimber(General):
    """
    The Hill Climber class optimizes a protein structure.
    It randomly changes a valid value in the protein structure.
    Each improvement or equivalent solution is kept for the next iteration.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).

    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    Notes
    -----
    This algorithm has built-in support for Simulated Annealing.

    Example
    -------
    >>> hill_climber = HillClimber("HHPHHHPH", 2)
    >>> hill_climber.run(show_plot = True, repeats = 5, iterations = 1000)

    """

    def __init__(self, protein_sequence: str, dimension: int) -> None:
        super().__init__(protein_sequence, dimension)

    def run(
        self,
        show_plot: bool = False,
        save_plot: bool = False,
        save_data: bool = False,
        repeats: int = 1,
        iterations: int = 1000
        ) -> None:
        """
        Runs the Hill Climber algorithm to optimize the protein structure.

        Parameters
        ----------
        show_plot : bool, optional
            If `True` show the plot. Default is `False`.

        save_plot : bool, optional
            If `True` save the plot. Default is `False`.

        save_data : bool, optional
            If `True`, saves the optimization results to a file.
            Default is `False`.

        repeats : int, optional
            The number of independent runs to perform. Default is `1`.

        iterations : int, optional
            The number of iterations per run. Default is `1000`.
        """
        self.run_algorithm(
            algorithm = "Hill Climber",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._hill_climber,
            accept_function = self._accept_function)

    def _hill_climber(
            self,
            iterations: int,
            accept_function: Callable[[int, int, float], tuple[bool, float]],
            temperature: float
            ) -> None:
        """Main function for running the Hill Climber algorithm."""
        best_score_list: list[int] = []
        for _ in range(iterations):
            best_score, protein, score_progression_list = self._hill_climber_helper(
                                                                    temperature,
                                                                    accept_function)
            # Save best results
            if best_score < self.best_score:
                self.best_protein = protein
                self.score_progression_list = score_progression_list
                self.best_score = best_score
            best_score_list.append(best_score)
            
        self.histogram_data.append(best_score_list)

    def _hill_climber_helper(
        self,
        temperature: float,
        accept_solution: Callable[[int, int, float], tuple[bool, float]],
        ) -> tuple[int, Protein, list[int]]:
        """
        Helper function that iteratively changes the protein structure
        and evaluated the results.
        It accepts a change if the score of the protein
        improved or stayed the same.
        The iterations stop when the score has not changed over 600 iterations.
        It has optional arguments for Simulated Annealing.
        It returns the best score,
        the `Protein` object and the list for scores
        for a score progression plot.

        Notes
        -----
        There is a restriction on the first protein structure generated.
        It forces to be a valid sequence.
        Without this restriction most structures would still be invalid,
        even after the algorithm.
        """
        score_progression_list: list[int] = []
        same_score_index: int = 0
        amino_directions: list[int] = random_fold(self.protein_sequence, self.dimension)
        protein: Protein = Protein(self.protein_sequence, amino_directions)
        protein.build_no_function()

        # Force valid solution
        while protein.protein_rating == 1:
            amino_directions = random_fold(self.protein_sequence, self.dimension)
            protein = Protein(self.protein_sequence, amino_directions)
            protein.build_no_function()

        best_rating: int = protein.protein_rating

        while same_score_index < MAX_UNCHANGED_ITERATIONS:
            # Choose random amino acid in sequence and give it new direction
            index = rd.randrange(len(amino_directions) - 1)

            if self.dimension == 2:
                new_direction = rd.choice([direction
                                           for direction in DIRECTION_CHOICES_2D
                                           if direction != amino_directions[index]])

            elif self.dimension == 3:
                new_direction = rd.choice([direction
                                           for direction in DIRECTION_CHOICES_3D
                                           if direction != amino_directions[index]])

            # Copy amino directions with new direction added
            new_amino_directions: list[int] = amino_directions[:]
            new_amino_directions[index] = new_direction

            # Create new Protein object and sequence
            protein_candidate: Protein = Protein(self.protein_sequence, new_amino_directions)
            protein_candidate.build_no_function()
            candidate_score = protein_candidate.protein_rating

            # Check for changes in score
            if candidate_score >= best_rating:
                same_score_index += 1
            elif candidate_score < best_rating:
                same_score_index = 0
            
            # Accept solution functions
            (accepted, temperature) = accept_solution(candidate_score, best_rating, temperature)
            if accepted:
                best_rating = candidate_score
                protein = protein_candidate
                amino_directions = new_amino_directions

            score_progression_list.append(best_rating)
        return best_rating, protein, score_progression_list
    
    def _accept_function(
            self,
            candidate_score: int,
            best_rating: int,
            temperature: float
            ) -> tuple[bool, float]:
        """
        Helper function for accepting the new rating.
        Returns `True` if smaller or equal.
        """
        if candidate_score <= best_rating:
            return True, temperature
        return False, temperature