import random as rd
import csv
import src.visualisation.plot_functions as plot

from src.classes.protein import Protein
from src.brute_force.brute_force import _direction_translator

class Random():
    """The Random class generates a random sequence for the folding direction"""

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.score_list: list[int] = []
        self.best_score: int = 0

    def run(self, show_plot: bool = False, save_plot: bool = False, n: int = 10000, save_data: bool = False) -> None:
        """Randomly generates sequences for a protein and calculates the scores."""
        self._random_iterated(show_plot, save_plot, n)
        if save_data:
            self.output_csv()

    def _random_iterated(self, show_plot: bool, save_plot: bool, n: int) -> None:
        """
        Iterates over multiple random generated folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.
        """
        for _ in range(n):
            protein = Protein(self.protein_sequence)
            protein.build_structure(self._random_fold)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein.build_structure(self._random_fold)

            self.score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score = protein.protein_rating
                best_protein = protein

        protein_sequence_map = {"HHPHHHPHPHHHPH" : "1",
                                "HPHPPHHPHPPHPHHPPHPH" : "2",
                                "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
                                "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
                                "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
                                "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
                                "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
                                "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

        folder = protein_sequence_map[best_protein.protein_sequence]
        base_path = "data/protein_random_folds/"

        # Plot histogram and visualize protein
        plot.histogram(self.protein_sequence, self.score_list, n, show = show_plot, save = save_plot, file_path = f"{base_path}{folder}", algorithm = "Random")
        plot.visualize(best_protein, show = show_plot, save = save_plot, file_path = f"{base_path}{folder}/best_random_fold")
        best_protein.output_csv(f"{base_path}{folder}/output")

    def _random_fold(self, protein_sequence: str) -> list[int]:
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
    
    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/random_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow(self.score_list)

        csvfile.close()