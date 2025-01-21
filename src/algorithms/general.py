from src.classes import Protein
from src.utils.helpers import save_and_visualize_results
from typing import Optional

class General():
    """
    A general template for protein structure optimization algorithms.
    The `General` class contains most basic functions every algorithm class would need.
    Every class that inherits `General` should implement their own `algorithm function`
    and then call `run_algorithm`.

    Parameter
    ---------
    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).
    
    Notes
    -----
    This class is not intended to be used directly.
    """

    def __init__(self, protein_sequence: str) -> None:
        """
        Attributes
        ----------
        `protein_squence` : str
            Protein sequence (for example `HHHPPPHPCCP`).

        `histogram_data` : list[list[`int`]]
            List containing all the scores for multiple runs.

        `best_score` : int
            The best score found. Default set to `0`.
        
        `best_protein` : `Protein`, optional
            The best `Protein` object found. Default set to `None`.
        
        `score_progression_list` : list[`int`]
            The score progression over iterations.
        """
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
        algorithm_function,
        check_solution_function = None,
        temperature: float = 2
        ) -> None:
        """
        Runs a protein optimization algorithm.

        Parameters
        ----------
        algorithm : str
            The name of the algorithm used (for example `Hill Climber`).
        
        show_plot : bool, optional
            If `True` show the plot. Default is `False`.

        save_plot : bool, optional
            If `True` save the plot. Default is `False`.

        save_data : bool, optional
            If `True`, saves the optimization results to a file. Default is `False`.

        repeats : int, optional
            The number of independent runs to perform. Default is `1`.

        iterations : int, optional
            The number of iterations per run. Default is `1000`.
        
        algortihm_function : Callable
            The main function for implementing the optimization algorithm.
            It updates the attributes initialized and has optional returns 
            of the form of a tuple consisting of best score, best protein and score progression list.
        
        check_solution_function: Callable, optional
            A function for that determines whether a new solution is accepted based on an acceptance probability.
            Default is None.
        
        temperature : int
            The initial temperature for annealing algorithm. Default is `2`.
        
        """
        for _ in range(repeats):
            if algorithm in ["Hill Climber", "Simulated Annealing"]:
                best_score_list: list[int] = []
                for _ in range(iterations):
                    best_score, protein, score_progression_list = algorithm_function(temperature, check_solution_function)
                    # Save best results
                    if best_score < self.best_score:
                        self.best_protein = protein
                        self.score_progression_list = score_progression_list
                        self.best_score = best_score
                    best_score_list.append(best_score)
                
                self.histogram_data.append(best_score_list)
            else:
                algorithm_function(iterations)
        
        # Save and visualize protein
        if self.best_protein is not None:
            save_and_visualize_results(
                best_protein = self.best_protein,
                algorithm = algorithm,
                histogram_data = self.histogram_data, 
                histogram = self.histogram_data[-1],
                iterations = iterations,
                show_plot= show_plot,
                save_plot= save_plot,
                save_data= save_data,
                score_progression = self.score_progression_list)
        else:
            print("Error: Did not find a valid protein.")
