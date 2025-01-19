import random as rd
import numpy as np

from src.utils.helpers import save_and_visualize_results
from .hill_climber import HillClimber

class SimulatedAnnealing(HillClimber):
    """
    The Simulated Annealing class optimizes a solution by mimicking the process of annealing in metals.
    It explores the solution space by accepting worse solutions with decreasing probability as the iterations progress.
    The algorithm balances exploration and exploitation to find an optimal solution.
    """

    def __init__(self, protein_sequence: str, temperature: int = 8):
        # Use init from Hill Climber class
        super().__init__(protein_sequence)

        # Initiate current temperature
        self.T: float = temperature
    
    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 1000) -> None:
        """Use hill climber algorithm with temperature to improve the sequence."""
        for _ in range(repeats):
            best_score_list: list[int] = []
            for _ in range(iterations):
                best_score, protein, score_progression_list = self._hill_climber(self._check_solution)
                # Save best results
                if best_score < self.best_score:
                    self.best_protein = protein
                    self.score_progression_list = score_progression_list
                    self.best_score = best_score
                best_score_list.append(best_score)
            
            self.histogram_data.append(best_score_list)

        # Save and visualize protein
        if self.best_protein is not None:
            save_and_visualize_results(self.best_protein, algorithm = "Simulated Annealing", histogram_data = self.histogram_data, 
            histogram = self.histogram_data[-1], iterations = iterations, show_plot= show_plot, save_plot= save_plot, save_data= save_data, score_progression = self.score_progression_list)
        else:
            print("Error: Did not find a valid protein.")

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