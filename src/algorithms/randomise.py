import random as rd
import src.visualisation.plot_functions as plot

from src.classes.protein import Protein
from src.brute_force.brute_force import _direction_translator
from typing import Optional

class Random():
    """The Random class generates a random sequence for the folding direction"""

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.best_score: int = 0
        self.best_protein: Optional[Protein] = None

    def run(self, show_plot: bool = False, save_plot: bool = False, n: int = 10000):
        self._random_iterated(show_plot, save_plot, n)

    def _random_iterated(self, show_plot: bool, save_plot: bool, n: int) -> None:
        """
        Iterates over multiple random generated folding sequences for a given protein string.
        Plots the distribution of the scores in a histogram.
        Shows the best structure that the random algorithm has generated.
        """
        score_list: list[int] = []

        for _ in range(n):
            protein = Protein(self.protein_sequence, self._random_fold)
            # Constraint to make sure a valid sequence is created
            while protein.protein_rating == 1:
                protein = Protein(self.protein_sequence, self._random_fold)

            score_list.append(protein.protein_rating)

            if protein.protein_rating < self.best_score:
                self.best_score = protein.protein_rating
                self.best_protein = protein

        protein_sequence_map = {"HHPHHHPHPHHHPH" : "1",
                                "HPHPPHHPHPPHPHHPPHPH" : "2",
                                "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
                                "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
                                "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
                                "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
                                "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
                                "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

        folder = protein_sequence_map[self.best_protein.protein_sequence]
        plot.histogram(protein, score_list, n, show = show_plot, save = save_plot, file_path = f"data/protein_random_folds/{folder}")

        if self.best_protein is not None:
            plot.visualize(self.best_protein, show = show_plot, save = save_plot, file_path = f"data/protein_random_folds/{folder}/best_random_fold") 

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