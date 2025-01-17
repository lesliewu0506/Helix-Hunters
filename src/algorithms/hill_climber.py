import random as rd
import csv
import copy
import src.visualisation.plot_functions as plot

from src.classes.protein import Protein
from src.algorithms.randomise import random_fold
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
        self.best_score_list: list[int] = []

        self.protein_sequence_map = {"HHPHHHPHPHHHPH" : "1",
                                    "HPHPPHHPHPPHPHHPPHPH" : "2",
                                    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
                                    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
                                    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
                                    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
                                    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
                                    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}
        self.folder = self.protein_sequence_map[protein_sequence]

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 1000) -> None:
        """Uses hill climbing algorithm to improve a random generated sequence."""
        for _ in range(repeats):
            best_score_list: list[int] = []
            for _ in range(iterations):
                best_score, protein, score_list = self._hill_climber()
                # Save best results
                if best_score < self.best_score:
                    self.best_protein = protein
                    self.best_score_list = score_list
                    self.best_score = best_score
                best_score_list.append(best_score)
            
            self.histogram_data.append(best_score_list)

        # Plot and save best protein structure
        base_path = "data/protein_hill_folds/"

        if self.best_protein is not None:
            plot.hill_visualizer(self.protein_sequence, self.best_score_list, show_plot = show_plot, save_plot = save_plot, file_path = f"{base_path}{self.folder}")
            plot.histogram(self.protein_sequence, self.histogram_data[-1], iterations = iterations, show = show_plot, save = save_plot, file_path = f"data/histogram_data/{self.folder}", algorithm = "Hill Climber")
            plot.visualize(self.best_protein, show = show_plot, save = save_plot, file_path = f"{base_path}{self.folder}/best_hill_fold")
            self.best_protein.output_csv(f"{base_path}{self.folder}/output")

        if save_data:
            self.output_csv()
    
    def _hill_climber(self, check_solution: Optional[Callable[[int, int], bool]] = None) -> tuple[int, Protein, list[int]]:
        score_list: list[int] = []
        same_score_index: int = 0
        amino_directions: list[int] = random_fold(self.protein_sequence)
        protein: Protein = Protein(self.protein_sequence, amino_directions)
        protein.build_no_function()

        # Force start with valid sequence
        while protein.protein_rating == 1:
            amino_directions = random_fold(self.protein_sequence)
            protein = Protein(self.protein_sequence, amino_directions)
            protein.build_no_function()

        best_rating: int = protein.protein_rating

        while same_score_index < 300:
            # Choose random amino acid in sequence and give it new direction
            index = rd.randrange(len(amino_directions))
            new_direction = rd.choice([direction for direction in [-2, -1, 1, 2] if direction != amino_directions[index]])

            # Copy amino directions with new direction added
            new_amino_directions: list[int] = copy.deepcopy(amino_directions)
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

            score_list.append(best_rating)
        return best_rating, protein, score_list

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/hill_climber_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        csvfile.close()