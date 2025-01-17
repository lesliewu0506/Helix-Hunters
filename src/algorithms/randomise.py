import random as rd
import csv
import src.visualisation.plot_functions as plot

from src.classes.protein import Protein
from src.brute_force.Brute_Force import _direction_translator

class Random():
    """The Random class generates a random sequence for the folding direction"""

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: list[list[int]] = []
        self.best_score: int = 0
        self.best_protein: Protein | None = None

        self.protein_sequence_map = {"HHPHHHPHPHHHPH" : "1",
                                    "HPHPPHHPHPPHPHHPPHPH" : "2",
                                    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
                                    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
                                    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
                                    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
                                    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
                                    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}
        self.folder = self.protein_sequence_map[protein_sequence]

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """Randomly generates sequences for a protein and calculates the scores."""
        for _ in range(repeats):
            self._random_iterated(iterations)

        base_path = "data/protein_random_folds/"
        # Plot histogram and visualize protein
        if self.best_protein is not None:
            plot.histogram(self.protein_sequence, self.histogram_data[-1], iterations, show = show_plot, save = save_plot, file_path = f"{base_path}{self.folder}", algorithm = "Random")
            plot.visualize(self.best_protein, show = show_plot, save = save_plot, file_path = f"{base_path}{self.folder}/best_random_fold")
            self.best_protein.output_csv(f"{base_path}{self.folder}/output")

        if save_data:
            self.output_csv()

    def _random_iterated(self, n: int) -> None:
        """
        Iterates over multiple random generated folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(random_fold)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein.build_structure(random_fold)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score = protein.protein_rating
                self.best_protein = protein

        self.histogram_data.append(score_list)

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/random_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        csvfile.close()

def random_fold(protein_sequence: str) -> list[int]:
    """
    Generates a random folding sequence.
    Uses relative directions (0, 1, 2) and translates them 
    into absolute directions (-2, -1, 1, 2).
    Returns the sequence as a list.
    """
    relative_direction_list: list[int] = []
    random_choice = [0, 1, 2]

    for _ in range(len(protein_sequence) - 2):
        direction = rd.choice(random_choice)
        relative_direction_list.append(direction)

    return _direction_translator(relative_direction_list)