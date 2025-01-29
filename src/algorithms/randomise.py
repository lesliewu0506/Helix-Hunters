from . import General
from src.classes import Protein
from src.utils import random_fold

class Random(General):
    """
    The Random class generates random sequences for a protein structure.
    It will run for `repeats` times with `iterations` iterations.
    At the end it will return the best protein structure found.
    Has optional arguments for showing and saving the data created by the run.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).

    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    Example
    -------
    >>> random = Random("HHPHHHPH", 2)
    >>> random.run(show_plot = True, repeats = 5, iterations = 1000)
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
        Runs the Random algorithm to generate protein structures.

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
            algorithm = "Random",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._random_iterated)

    def _random_iterated(
            self,
            n: int,
            accept_function: None,
            temperature: float
            ) -> None:
        """
        Helper function that generates multiple random
        folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.

        Notes
        -----
        - Only valid structures are accepted.
        Without this restriction,
        most of the generated structures would be invalid.
        - `accept_function` and `temperature` are dummy variables,
        that are not intended to be used. 
        They are there to be compatible with the General class.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(random_fold, self.dimension)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein.build_structure(random_fold, self.dimension)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score: int = protein.protein_rating
                self.best_protein = protein

        self.histogram_data.append(score_list)