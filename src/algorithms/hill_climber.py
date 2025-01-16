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
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Protein | None = None
        self.score_list: list[int] = []

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
            self._hill_climber(show_plot, save_plot, iterations)
        
        if save_data:
            self.output_csv()
    
    def _hill_climber(self, show_plot: bool, save_plot: bool, n: int):
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
            self.score_list.append(best_rating)

        plt.plot(list(range(1, n + 1)), self.score_list)
        plt.show()

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/hill_climber_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        csvfile.close()