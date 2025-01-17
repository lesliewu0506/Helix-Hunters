import random as rd
import csv
import numpy as np
import src.visualisation.plot_functions as plot

from .hill_climber import HillClimber

class SimulatedAnnealing(HillClimber):
    """
    The Simulated Annealing class optimizes a solution by mimicking the process of annealing in metals.
    It explores the solution space by accepting worse solutions with decreasing probability as the iterations progress.
    The algorithm balances exploration and exploitation to find an optimal solution.
    """

    def __init__(self, protein_sequence: str, temperature: int = 3):
        # Use init from Hill Climber class
        super().__init__(protein_sequence)

        # Initiate current temperature
        self.T: float = temperature
    
    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 1000) -> None:
        """Use hill climber algorithm with temperature to improve the sequence."""
        for _ in range(repeats):
            best_score_list: list[int] = []
            for _ in range(iterations):
                best_score, protein, score_list = self._hill_climber(self.check_solution)
                # Save best results
                if best_score < self.best_score:
                    self.best_protein = protein
                    self.best_score_list = score_list
                    self.best_score = best_score
                best_score_list.append(best_score)
            
            self.histogram_data.append(best_score_list)
        # Plot and save best protein structure
        base_path = "data/protein_annealing_folds/"

        if self.best_protein is not None:
            plot.hill_visualizer(self.protein_sequence, self.best_score_list, show_plot = show_plot, save_plot = save_plot, file_path = f"{base_path}{self.folder}", algorithm = "Simulated Annealing")
            plot.histogram(self.protein_sequence, self.histogram_data[-1], iterations = iterations, show = show_plot, save = save_plot, file_path = f"{base_path}{self.folder}", algorithm = "Simulated Annealing")
            plot.visualize(self.best_protein, show = show_plot, save = save_plot, file_path = f"{base_path}{self.folder}/best_annealing_fold")
            self.best_protein.output_csv(f"{base_path}{self.folder}/output")

        if save_data:
            self.output_csv()

    def check_solution(self, new_rating: int, old_rating: int) -> bool:
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
        alpha: float = 0.999
        self.T = self.T * alpha

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/annealing_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        csvfile.close()