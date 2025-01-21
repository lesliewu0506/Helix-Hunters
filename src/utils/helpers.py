"""
This module contains general-purpose utility functions that can be used across 
the project. These functions are designed to simplify repetitive tasks, improve 
code readability, and promote reusability.

Helper functions in this file are not specific to any single module but can be 
used by various parts of the project wherever needed.
"""
import random as rd
import csv

import src.visualisation.plot_functions as plot
from src.classes import Protein
from src.utils.constants import protein_sequence_map, algorithm_folder_map

# ===============================================================
# Utility Functions
# ===============================================================

def random_fold(protein_sequence: str) -> list[int]:
    """
    Generates a random folding sequence using relative direction (0, 1, 2).
    Translates them into absolute direction (-2, -1, 1, 2).

    Parameter
    ----------
    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).

    Returns
    -------
    list[`int`]
        A list of absolute directions.

    Notes
    -----
    This function does not check if the sequence collides with itself.
    A seperate function is needed to check this.
    """
    relative_direction_list: list[int] = []
    random_choice = [0, 1, 2]

    for _ in range(len(protein_sequence) - 2):
        direction = rd.choice(random_choice)
        relative_direction_list.append(direction)

    return direction_translator(relative_direction_list)

def direction_translator(directions: list[int]) -> list[int]:
    """
    Translates relative paths (0, 1, 2) to absolute paths (-2, -1, 1, 2)
    Also adds 1 and 0 at start and end respectively for correct format.
    
    Parameter
    ---------
    directions : list[`int`]
        A list of relative directions.

    Returns
    -------
    list[`int`]
        A list of absolute directions.
    """
    direction_map: dict[int, list[int]] = {1: [2, 1, -2], -1: [-2, -1, 2], 2: [-1, 2, 1], -2: [1, -2, -1]}
    folding_sequence: list[int] = [1]

    for i, direction in enumerate(directions):
        folding_sequence.append(direction_map[folding_sequence[i]][direction])

    folding_sequence.append(0)
    return folding_sequence

# ===============================================================
# Visualize Helpers
# ===============================================================

def save_and_visualize_results(
    best_protein: Protein,
    algorithm: str,
    histogram_data: list[list[int]], 
    histogram: list[int], 
    iterations: int, 
    show_plot: bool, 
    save_plot: bool, 
    save_data: bool, 
    score_progression: list[int] = []
    ) -> None:
    """
    Main function for visualizing and saving the results of the optimization algorithm.
    It plots the histogram for the distribution of scores of an algorithm for a protein sequence.
    It also plots the best protein structure found from that particular algorithm.
    Its structure is then written into a `CSV` file with a specific format.
    If the algorithm is either `Hill Climber` or `Simulated Annealing`, it will also plot
    and save the score progression plot.
    It can also write the distribution of scores to a `CSV` file.

    Parameters
    ----------
    best_protein : `Protein` object
        Best protein object containing the sequence and folding structure.

    algorithm : str
        The name of the algorithm used (for example `Hill Climber`).

    histogram_data : list[list[`int`]]
        List containing all the scores for multiple runs.

    histogram : list[int]
        Histogram data for last run.

    iterations : int
        Amount of iterations for the run.

    show_plot : bool
        Whether to show the plots.

    save_plot : bool
        Whether to save the plots.

    save_data : bool
        Whether to save histogram data to a `CSV` file.

    score_progression: list[`int`], optional
        The score progression over iterations. Default is empty list `[]`.
    """
    protein_sequence = best_protein.protein_sequence

    # Create correct path for each algorithm and sequence
    folder = protein_sequence_map[protein_sequence]
    base_path: str = f"data/protein_{algorithm_folder_map[algorithm]}_folds/{folder}"

    # Plots the progression of Hill Climber/Simulated Annealing algorithm
    if algorithm in ["Hill Climber", "Simulated Annealing"]:
        plot.score_progression(protein_sequence, score_progression, show = show_plot, save = save_plot, file_path = f"{base_path}", algorithm = algorithm)
    
    # Plots score distribution for one repeat of algorithm
    plot.histogram(protein_sequence, histogram, iterations = iterations, show = show_plot, save = save_plot, file_path = f"{base_path}", algorithm = algorithm)

    # Plots the best protein structure for one repeat of algorithm
    plot.visualize_protein(best_protein, algorithm, show = show_plot, save = save_plot, file_path = f"{base_path}")

    if save_data:
        output_histogram_csv(protein_sequence, algorithm, histogram_data)
        best_protein.output_csv(f"{base_path}/output")

# ===============================================================
# File Operations
# ===============================================================

def output_histogram_csv(
    protein_sequence: str,
    algorithm: str,
    histogram_data: list[list[int]]
    ) -> None:
    """
    Main function for saving the histogram data into a `CSV` file.
    It writes the data into a designated location.
    
    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).

    algorithm : str
        The name of the algorithm used (for example `Hill Climber`).

    histogram_data : list[list[`int`]]
        List containing all the scores for multiple runs.
    """
    folder = protein_sequence_map[protein_sequence]
    with open(f"data/histogram_data/{folder}/{algorithm.title()}_{protein_sequence}.csv", 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile)

        for histogram in histogram_data:
            writer.writerow(histogram)

    csvfile.close()