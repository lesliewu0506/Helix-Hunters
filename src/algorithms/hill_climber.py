import random as rd

from . import General
from src.utils import random_fold
from src.classes import Protein
from typing import Callable, Optional

class HillClimber(General):
    """
    The Hill Climber class generates a sequence for the folding direction.
    Randomly changes a valid value. 
    Each improvement or equivalent solution is kept for the next iteration.
    """

    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence)

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 1000) -> None:
        """Uses hill climbing algorithm to improve a random generated sequence."""
        self.run_algorithm(
            algorithm = "Hill Climber",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._hill_climber)
    
    def _hill_climber(self, temperature: float = 1, check_solution: Optional[Callable[[int, int, float], tuple[bool, float]]] = None) -> tuple[int, Protein, list[int]]:
        score_progression_list: list[int] = []
        same_score_index: int = 0
        amino_directions: list[int] = random_fold(self.protein_sequence)
        protein: Protein = Protein(self.protein_sequence, amino_directions)
        protein.build_no_function()

        # Force valid solution
        while protein.protein_rating == 1:
            amino_directions = random_fold(self.protein_sequence)
            protein = Protein(self.protein_sequence, amino_directions)
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
                (accepted, temperature) = check_solution(candidate_score, best_rating, temperature)
                if accepted:
                    best_rating = candidate_score
                    protein = protein_candidate
                    amino_directions = new_amino_directions

            score_progression_list.append(best_rating)
        return best_rating, protein, score_progression_list