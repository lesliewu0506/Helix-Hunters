from src.classes import Protein
from src.utils import save_and_visualize_results
from typing import Optional

class General():

    def __init__(self, protein_sequence: str) -> None:
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
        temperature: float = 1
        ) -> None:

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
