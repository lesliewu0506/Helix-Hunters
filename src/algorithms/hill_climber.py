import random as rd
import csv
import copy
import matplotlib.pyplot as plt
import src.visualisation.plot_functions as plot

from src.classes.protein import Protein
from src.algorithms.randomise import random_fold
class HillClimber():
    """
    The Hill Climber class generates a sequence for the folding direction.
    Randomly changes a valid value. 
    Each improvement or equivalent solution is kept for the next iteration.
    """

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[int] = []
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

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 3000) -> None:
        """Uses hill climbing algorithm to improve a random generated sequence."""
        for _ in range(repeats):
            best_score, protein, score_list = self._hill_climber(iterations)
            # Save best results
            if best_score < self.best_score:
                self.best_protein = protein
                self.best_score_list = score_list
                self.best_score = best_score

            self.histogram_data.append(best_score)

        # Plot and save best protein structure
        plot.hill_visualizer(self.protein_sequence, self.best_score_list, show_plot = show_plot, save_plot = save_plot, file_path = f"data/protein_hill_folds/{self.folder}")
        plot.histogram(self.protein_sequence, self.histogram_data, iterations = iterations, show = show_plot, save = save_plot, file_path = f"data/histogram_data/{self.folder}", algorithm = "hill_climber")
        plot.visualize(self.best_protein, show = show_plot, save = save_plot, file_path = f"data/protein_hill_folds/{self.folder}/best_hill_fold")

        if save_data:
            self.output_csv()
    
    def _hill_climber(self, n: int) -> tuple[int, Protein, list[int]]:
        score_list: list[int] = []

        amino_directions: list[int] = random_fold(self.protein_sequence)
        protein: Protein = Protein(self.protein_sequence, amino_directions)
        protein.build_no_function()

        # Force start with valid sequence
        while protein.protein_rating == 1:
            amino_directions: list[int] = random_fold(self.protein_sequence)
            protein: Protein = Protein(self.protein_sequence, amino_directions)
            protein.build_no_function()

        best_rating: int = protein.protein_rating

        for _ in range(n):
            # Choose random amino acid in sequence and give it new direction
            index = rd.randrange(len(amino_directions))
            new_direction = rd.choice([direction for direction in [-2, -1, 1, 2] if direction != amino_directions[index]])

            # Copy amino directions with new direction added
            new_amino_directions: list[int] = copy.deepcopy(amino_directions)
            new_amino_directions[index] = new_direction

            # Create new Protein object and sequence
            protein_candidate: Protein = Protein(self.protein_sequence, new_amino_directions)
            protein_candidate.build_no_function()

            # Accept candidate if score becomes lower
            if protein_candidate.protein_rating <= best_rating:
                best_rating = protein_candidate.protein_rating
                protein = protein_candidate
                amino_directions = new_amino_directions
            score_list.append(best_rating)

        return best_rating, protein, score_list

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/hill_climber_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(self.histogram_data)

        csvfile.close()