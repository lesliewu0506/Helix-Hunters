import random as rd

from src.utils.helpers import save_and_visualize_results, random_fold
from src.classes.protein import Protein
from typing import Callable, Optional

class HillClimber():
    """
    The Hill Climber class generates a sequence for the folding direction.
    Randomly changes a valid value. 
    Each improvement or equivalent solution is kept for the next iteration.
    """

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Protein | None = None
        self.score_progression_list: list[int] = []

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """Uses hill climbing algorithm to improve a random generated sequence."""
        for _ in range(repeats):
            best_score_list: list[int] = []
            for _ in range(iterations):
                best_score, protein, score_progression_list = self._hill_climber()
                # Save best results
                if best_score < self.best_score:
                    self.best_protein = protein
                    self.score_progression_list = score_progression_list
                    self.best_score = best_score
                best_score_list.append(best_score)
            
            self.histogram_data.append(best_score_list)

        # Save and visualize protein
        if self.best_protein is not None:
            save_and_visualize_results(self.best_protein, algorithm = "Hill Climber", histogram_data = self.histogram_data, 
            histogram = self.histogram_data[-1], iterations = iterations, show_plot= show_plot, save_plot= save_plot, save_data= save_data, score_progression = self.score_progression_list)
        else:
            print("Error: Did not find a valid protein.")
    
    def _hill_climber(self, check_solution: Optional[Callable[[int, int], bool]] = None) -> tuple[int, Protein, list[int]]:
        score_progression_list: list[int] = []
        same_score_index: int = 0
        amino_directions: list[int] = random_fold(self.protein_sequence)
        protein: Protein = Protein(self.protein_sequence, amino_directions)
        protein.build_no_function()

        # Force valid solution
        while protein.protein_rating == 1:
            amino_directions = random_fold(self.protein_sequence)
            protein: Protein = Protein(self.protein_sequence, amino_directions)
            protein.build_no_function()
        best_rating: int = protein.protein_rating

        while same_score_index < 600:
            # Choose random amino acid in sequence and give it new direction
            index = rd.randrange(len(amino_directions))
            new_direction = rd.choice([direction for direction in [-2, -1, 1, 2] if direction != amino_directions[index]])

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
            
            # Check for simulated annealing
            if check_solution is None:
                # Accept candidate if score becomes lower
                if candidate_score <= best_rating:
                    best_rating = candidate_score
                    protein = protein_candidate
                    amino_directions = new_amino_directions
            else:
                # Accept based on probability
                if check_solution(candidate_score, best_rating):
                    best_rating = candidate_score
                    protein = protein_candidate
                    amino_directions = new_amino_directions

            score_progression_list.append(best_rating)
        return best_rating, protein, score_progression_list