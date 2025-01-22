import random as rd
import numpy as np

from . import HillClimber

class SimulatedAnnealing(HillClimber):
    """
    The Simulated Annealing class optimizes a protein structure by mimicking the process of annealing in metals.
    It searches for multiple solutions by accepting worse solutions with decreasing probability as the iterations progress.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).
    
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    temperature : int, optional
        The initial temperature for annealing algorithm. Default is `2`.

    Notes
    -----
    This class inherits most of the `HillClimber` class functions. 
    """

    def __init__(self, protein_sequence: str, dimension: int, temperature: int = 2) -> None:
        super().__init__(protein_sequence, dimension)

        self.T: float = temperature

    def run(
        self,
        show_plot: bool = False,
        save_plot: bool = False,
        save_data: bool = False,
        repeats: int = 1,
        iterations: int = 1000
        ) -> None:
        """
        Runs the Simulated Annealing algorithm to optimize the protein structure.

        Parameters
        ----------
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
        """
        self.run_algorithm(
            algorithm = "Simulated Annealing",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._hill_climber,
            check_solution_function = self._check_solution,
            temperature= self.T)

    def _check_solution(self, new_rating: int, old_rating: int, temperature: float) -> tuple[bool, float]:
        """
        Helper function that determines whether a new solution is accepted based on the acceptance probability.
        The acceptance probability is calculated using the Metropolis criterion, where 
        worse solutions are accepted with a probability that decreases as the temperature lowers.
        """
        delta: int = new_rating - old_rating

        if delta <= 0:
            # Update temperature
            temperature = self._update_temperature(temperature)
            return True, temperature
        
        probability: float = np.exp(-delta / temperature)

        # Update temperature
        temperature = self._update_temperature(temperature)
        # Return if accepted new rating
        if rd.random() < probability:
            return True, temperature
        else: 
            return False, temperature

    def _update_temperature(self, temperature) -> float:
        """
        Updates the current temperature using an exponential decay formula.

        Notes
        -----
        This decay function drops extremely fast. This is due to the fact that
        higher values (0.7, 0.8, 0.9) causes many invalid protein structures.
        By trial and error, the value of 0.6 has been found for this project.
        """
        temperature = temperature * 0.6
        return temperature