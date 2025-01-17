"""
This module contains general-purpose utility functions that can be used across 
the project. These functions are designed to simplify repetitive tasks, improve 
code readability, and promote reusability.

Helper functions in this file are not specific to any single module but can be 
used by various parts of the project wherever needed.
"""
import csv
import src.visualisation.plot_functions as plot
from src.classes.protein import Protein

protein_sequence_map: dict[str, str] = {
    "HHPHHHPHPHHHPH" : "1",
    "HPHPPHHPHPPHPHHPPHPH" : "2",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

algorithm_folder_map: dict[str, str] = {
    "Random" : "random",
    "Greedy" : "greedy",
    "Hill Climber" : "hill",
    "Simulated Annealing" : "annealing",
    "Plant Propagation" : "propagation",
    "Genetic Algorithm" : "genetic"}

def save_and_visualize_results(best_protein: Protein, algorithm: str, histogram_data: list[list[int]], histogram: list[int], iterations: int, show_plot: bool, 
                               save_plot: bool, save_data: bool, score_progression: list[int] = []):
    """Saves and visualizes the results of the optimization algorithm."""
    protein_sequence = best_protein.protein_sequence

    # Create correct path for each algorithm and sequence
    folder = protein_sequence_map[protein_sequence]
    base_path: str = f"data/protein_{algorithm_folder_map[algorithm]}_folds/{folder}"

    # Plots the progression of Hill Climber/Simulated Annealing algorithm
    if algorithm in ["Hill Visualizer", "Simulated Annealing"]:
        plot.hill_visualizer(protein_sequence, score_progression, show_plot = show_plot, save_plot = save_plot, file_path = f"{base_path}", algorithm = algorithm)
    
    # Plots score distribution for one repeat of algorithm
    plot.histogram(protein_sequence, histogram, iterations = iterations, show = show_plot, save = save_plot, file_path = f"{base_path}", algorithm = algorithm)

    # Plots the best protein structure for one repeat of algorithm
    plot.visualize(best_protein, algorithm, show = show_plot, save = save_plot, file_path = f"{base_path}")

    if save_data:
        output_histogram_csv(protein_sequence, algorithm, histogram_data)
        best_protein.output_csv(f"{base_path}/output")

def output_histogram_csv(protein_sequence: str, algorithm: str, histogram_data: list[list[int]]) -> None:
    """Saves histogram data into a csv file."""
    folder = protein_sequence_map[protein_sequence]
    with open(f"data/histogram_data/{folder}/{algorithm}_{protein_sequence}.csv", 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile)

        for histogram in histogram_data:
            writer.writerow(histogram)

    csvfile.close()

def direction_translator(directions: list[int]) -> list[int]:
    """
    Helper function for translating the relative paths (0, 1, 2),
    to absolute paths (-2, -1, 1, 2). Returns list of directions.
    """
    direction_map: dict[int, list[int]] = {1: [2, 1, -2], -1: [-2, -1, 2], 2: [-1, 2, 1], -2: [1, -2, -1]}
    folding_sequence: list[int] = [1]

    for i, direction in enumerate(directions):
        folding_sequence.append(direction_map[folding_sequence[i]][direction])

    folding_sequence.append(0)
    return folding_sequence