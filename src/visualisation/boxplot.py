import csv
import matplotlib.pyplot as plt
import numpy as np

from src.utils.constants import algorithms, protein_sequence_map

def boxplot(protein_sequence: str) -> None:
    """
    Main function for plotting a boxplot of the score distributions
    of different algorithms of a protein sequence.
    It imports the data from CSV files that corresponds to the protein sequence.
    Then it couples each data with the corresponding algorithm in a `dict`. 
    The minimum score is then calculated for the limits for the boxplots.
    The plots are automatically saved in a designated folder.

    Parameters
    ----------
    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).
    """
    folder: str = protein_sequence_map[protein_sequence]
    data_structure: dict[str, list[int]] = _create_data_structure(protein_sequence, folder)
    data_list: list[list[int]] = [data for data in data_structure.values()]
    min_score: int = _get_minimum_score(data_structure)

    fig, ax = plt.subplots(figsize = (12, 7))
    box = ax.boxplot(data_list, patch_artist = True, labels = algorithms)

    # Customize boxplot appearance
    colors = plt.cm.viridis(np.linspace(0, 1, len(algorithms)))
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.85)

    ax.set_title(f"Protein Score Distribution for Different Algorithms\nProtein Sequence: {protein_sequence}", fontsize = 14, fontweight = "bold")
    ax.set_xlabel("Algorithm", fontsize = 12)
    ax.set_ylabel("Score", fontsize = 12)

    ax.set_ylim(min_score - 1 , 2)
    ax.grid(axis = "y", linestyle = "--", alpha = 0.85)

    plt.tight_layout()
    plt.savefig(f"data/results/Distributions For {protein_sequence}.png", dpi = 600)
    plt.show()

def _create_data_structure(protein_sequence: str, folder: str) -> dict[str, list[int]]:
    """
    Helper function that creates a dictionary
    with algorithm names as keys and the boxplot data as values.
    It returns the dictionary.
    """
    data: dict[str, list[int]] = {}
    for algorithm in algorithms:
        data[algorithm] = _import_data(protein_sequence, algorithm, folder)
    return data

def _import_data(protein_sequence: str, algorithm: str, folder: str) -> list[int]:
    """
    Helper function that imports the scores from a CSV file.
    It returns a list of the scores.
    """
    histogram_data: list[int] = []

    with open(f"data/histogram_data/{folder}/{algorithm}_{protein_sequence}.csv", "r") as csvfile:

        reader = csv.reader(csvfile)

        [[histogram_data.append(int(value)) for value in row] for row in reader]

    csvfile.close()

    return histogram_data

def _get_minimum_score(data: dict[str, list[int]]) -> int:
    """
    Helper function that calculates the lowest score of a given data structure.
    It returns the lowest score.
    """
    minimum_list: list[int] = []
    for histogram_data in data.values():
        minimum_list.append(min(histogram_data))
    return min(minimum_list)