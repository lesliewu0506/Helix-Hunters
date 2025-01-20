import random as rd
import numpy as np

from .hill_climber import HillClimber

class SimulatedAnnealing(HillClimber):
    """
    The Simulated Annealing class optimizes a solution by mimicking the process of annealing in metals.
    It explores the solution space by accepting worse solutions with decreasing probability as the iterations progress.
    The algorithm balances exploration and exploitation to find an optimal solution.
    """

    def __init__(self, protein_sequence: str, temperature: int = 8) -> None:
        # Use init from General Class
        super().__init__(protein_sequence)

        # Initiate current temperature
        self.T: float = temperature
    
    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 1000) -> None:
        """Use hill climber algorithm with temperature to improve the sequence."""
        self.run_algorithm(
            algorithm = "Simulated Annealing",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._hill_climber,
            check_solution_function = self._check_solution)

    def _check_solution(self, new_rating: int, old_rating: int) -> bool:
        """
        Calculates the acceptance rate of a new change.
        Returns True if accepted, else False.
        """
        delta: int = new_rating - old_rating

        if delta <= 0:
            # Update temperature
            self._update_temperature()
            return True
        
        probability: float = np.exp(-delta / self.T)

        # Update temperature
        self._update_temperature()
        # Return if accepted new rating
        if rd.random() < probability:
            return True
        else: 
            return False

    def _update_temperature(self) -> None:
        """Updates temperature based on exponential decay."""
        alpha: float = 0.97
        self.T = self.T * alpha