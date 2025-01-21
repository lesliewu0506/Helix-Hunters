import csv
import matplotlib.pyplot as plt
import numpy as np

algorithms: list[str] = ["random", "greedy", "Hill Climber", "Simulated Annealing"]

def boxplot(protein_sequence: str, folder: str):
    data_structure: dict[str, list[int]] = create_data_structure(protein_sequence, folder)
    algorithm_list: list[str] = list(data_structure.keys())
    data_list: list[list[int]] = [data for data in data_structure.values()]
    min_score: int = get_min(data_structure)
    
    fig, ax = plt.subplots(figsize = (12, 7))
    box = ax.boxplot(data_list, patch_artist = True, labels = algorithm_list)

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

def import_data(protein_sequence: str, algorithm: str, folder: str) -> list[int]:
    histogram_data: list[int] = []

    with open(f"data/histogram_data/{folder}/{algorithm}_{protein_sequence}.csv", "r") as csvfile:

        reader = csv.reader(csvfile)

        [[histogram_data.append(int(value)) for value in row] for row in reader]

    csvfile.close()

    return histogram_data

def create_data_structure(protein_sequence: str, folder: str) -> dict[str, list[int]]:

    data: dict[str, list[int]] = {}
    for algorithm in algorithms:
        data[algorithm] = import_data(protein_sequence, algorithm, folder)
    return data

def get_min(data: dict[str, list[int]]) -> int:
    minimum_list: list[int] = []
    for histogram_data in data.values():
        minimum_list.append(min(histogram_data))
    return min(minimum_list)