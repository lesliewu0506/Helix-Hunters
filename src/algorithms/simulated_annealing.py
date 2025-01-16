import random as rd
import csv

import src.visualisation.plot_functions as plot

from src.classes.protein import Protein

class SimulatedAnnealing():
    """
    The Simulated Annealing class optimizes a solution by mimicking the process of annealing in metals.
    It explores the solution space by accepting worse solutions with decreasing probability as the iterations progress.
    The algorithm balances exploration and exploitation to find an optimal solution.
    """

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
    
    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        with open(f"data/histogram_data/{self.folder}/annealing_{self.protein_sequence}.csv", 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        csvfile.close()