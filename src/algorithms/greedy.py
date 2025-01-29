import random as rd

from . import General
from src.classes import Protein
from src.utils import direction_translator, DIRECTIONS_2D, DIRECTIONS_3D

class Greedy(General):
    """
    The Greedy random class generates a sequence for the folding direction.
    It uses a greedy algorithm every 5 iterations to try and improve the score.
    
    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).
    
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    Example
    -------
    >>> greedy = Greedy("HHPHHHPH", 2)
    >>> greedy.run(show_plot = True, repeats = 5, iterations = 1000)
    """

    def __init__(self, protein_sequence: str, dimension: int) -> None:
        super().__init__(protein_sequence, dimension)

    def run(
        self,
        show_plot: bool = False,
        save_plot: bool = False,
        save_data: bool = False,
        repeats: int = 1,
        iterations: int = 10000
        ) -> None:
        """
        Runs the Greedy algorithm to generate protein structures.

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
            The number of iterations per run. Default is `10000`.
        """
        self.run_algorithm(
            algorithm = "Greedy",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._greedy_iterated)

    def _greedy_iterated(
            self,
            n: int,
            accept_function: None,
            temperature: float
            ) -> None:
        """
        Helper function that generates multiple protein structures
        for a given protein sequence.
        For each iteration,
        it will use a mix of random
        and greedy algorithm to create the structure.
        It saves the best protein and the distributions of scores.

        Notes
        -----
        There is a restriction on the first protein structure generated.
        It forces to be a valid sequence.
        Without this restriction most structures would still be invalid,
        even after the algorithm.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(self._greedy_fold, self.dimension)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein = Protein(self.protein_sequence)
                protein.build_structure(self._greedy_fold, self.dimension)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score: int = protein.protein_rating
                self.best_protein = protein

        self.histogram_data.append(score_list)

    def _greedy_fold(self, protein_sequence: str, dimension: int) -> list[int]:
        """
        Helper function that generates a random folding sequence.
        Every 5 iterations it will make a greedy decision.
        Uses relative directions (0, 1, 2, 3, 4) and translates them
        into absolute directions (-3, -2, -1, 1, 2, 3).
        Returns the sequence as a list.

        Notes
        -----
        Does not implement greedy decision on every iteration.
        This would cause the function to check every possible direction.
        The total run time would be enormous. That is why every fifth
        iteration a greedy decision is made.
        """
        relative_direction_list: list[int] = []

        if dimension == 2:
            directions = DIRECTIONS_2D

        elif dimension == 3:
            directions = DIRECTIONS_3D
        
        for i in range(len(protein_sequence) - 2):
            
            if i % 5 == 0:
                direction = self._check_greedy(relative_direction_list, dimension)
                relative_direction_list.append(direction)

            else:
                direction = rd.choice(directions)
                relative_direction_list.append(direction)

        return direction_translator(relative_direction_list, dimension)

    def _check_greedy(self, direction_list: list[int], dimension: int) -> int:
        """
        Helper function that implements a greedy algorithm
        to find the best possible next direction.
        Evaluates all possible continuations and
        selecting lowest score directions randomly if there are ties.
        """
        # Evaluate all possible next directions and collect scores
        protein_scores = []
        best_directions = []
        minimum_score = 1

        for direction in (DIRECTIONS_3D if dimension == 3 else DIRECTIONS_2D):

            new_direction_list = direction_list + [direction]
            abs_direction = direction_translator(new_direction_list, dimension)
            abs_direction.remove(0)

            # Create a Protein object and compute its rating
            protein = Protein(self.protein_sequence[:len(new_direction_list)],
                              amino_directions=abs_direction)
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