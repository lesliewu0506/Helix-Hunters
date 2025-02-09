from src.classes import Protein
from src.utils import save_and_visualize_results, TEMPERATURE, DIMENSIONS
from typing import Optional, Callable

class General():
    """
    A general template for protein structure optimization algorithms.
    The `General` class contains most basic functions
    every algorithm class would need.
    Every class that inherits `General`
    should implement their own `algorithm function`
    and then call `run_algorithm`.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).
    
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    Notes
    -----
    This class is not intended to be used directly.

    Example
    -------
    ```python
    class HillClimber(General):
        def __init__(self, protein_sequence: str, dimension: int) -> None:
            super().__init__(protein_sequence, dimension)

        def run(self, ...):
            self.run_algorithm(...)
    ```
    """

    def __init__(self, protein_sequence: str, dimension: int) -> None:
        """
        Attributes
        ----------
        `dimension` : int
            The dimension in which the folding takes place (`2` or `3`).

        `protein_squence` : str
            Protein sequence (for example `HHPHHHPH`).

        `histogram_data` : `list[list[int]]`
            List containing all the scores for multiple runs.

        `best_score` : int
            The best score found. Default set to `0`.
        
        `best_protein` : `Protein`, optional
            The best `Protein` object found. Default set to `None`.
        
        `score_progression_list` : `list[int]`
            The score progression over iterations.

        Raises
        ------
        ValueError
            If the dimension is not `2` or `3`.
        """
        if dimension not in DIMENSIONS:
            raise ValueError(f"Invalid dimension given. Choose from:\n{DIMENSIONS}.")

        self.dimension: int = dimension
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Optional[Protein] = None
        self.score_progression_list: list[int] = []

    def run_algorithm(
        self,
        algorithm: str,
        show_plot: bool,
        save_plot: bool,
        save_data: bool,
        repeats: int,
        iterations: int,
        algorithm_function: Callable[[int, Optional[Callable[[int, int, float], tuple[bool, float]]], float], None],
        accept_function: Optional[Callable[[int, int, float], tuple[bool, float]]] = None,
        temperature: float = TEMPERATURE
        ) -> None:
        """
        Runs a protein optimization algorithm.
        Can show and save the data created with arguments.

        Parameters
        ----------
        algorithm : str
            The name of the algorithm used (for example `Hill Climber`).
        
        show_plot : bool
            If `True` show the plot.

        save_plot : bool
            If `True` save the plot.

        save_data : bool
            If `True`, saves the optimization results to a file.

        repeats : int
            The number of independent runs to perform.

        iterations : int
            The number of iterations per run.
        
        algortihm_function : Callable
            The main function for implementing the optimization algorithm.
            It updates the attributes initialized.
        
        check_solution_function: Callable, optional
            A function for that determines whether a new solution
            is accepted based on an acceptance probability.
            Default is `None`.
        
        temperature : int, optional
            The initial temperature for annealing algorithm.
            Default is `TEMPERATURE`.

        Notes
        -----
        The method saves the best protein structure,
        score progression, and histogram data.
        These results can be visualized
        or saved to files for further analysis.

        Raises
        ------
        ValueError
            If repeats or iterations has an invalid value (<1).

        RunTimeError
            If no valid protein structure was found
            during the optimization process.
        """
        if repeats < 1 or iterations < 1:
            raise ValueError("Both repeats and iterations must be at least 1.")

        for _ in range(repeats):
            algorithm_function(
                iterations,
                accept_function = accept_function,
                temperature = temperature)

        # Save and visualize protein
        if self.best_protein is not None:
            save_and_visualize_results(
                dimension = self.dimension,
                best_protein = self.best_protein,
                algorithm = algorithm,
                histogram_data = self.histogram_data,
                histogram = self.histogram_data[-1],
                iterations = iterations,
                show_plot = show_plot,
                save_plot = save_plot,
                save_data = save_data,
                score_progression = self.score_progression_list)
        else:
            raise RuntimeError("No valid protein structure found during the optimization process.")