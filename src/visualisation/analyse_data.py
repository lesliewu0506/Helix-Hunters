import csv
import matplotlib.pyplot as plt
import numpy as np

algorithms: list[str] = ["random", "greedy", "Hill Climber", "Simulated Annealing"]

def histogram_3d(protein_sequence: str, folder: str):
    data: dict[str, list[list[int]]] = create_data_structure(protein_sequence, folder)
    min_score: int = get_min(data)

    bins = np.arange(min_score - 0.5, 1.5, 1)

    frequency_data: list[int] = []
    for algorithm in algorithms:
        # weights = np.ones_like(data[algorithm]) / len(data[algorithm])
        histogram, _ = np.histogram(data[algorithm], bins = bins, density = True)
        frequency_data.append(histogram)

    frequency_data = np.array(frequency_data)

    fig = plt.figure(figsize = (12, 7))
    ax = fig.add_subplot(111, projection = "3d")

    xpos, ypos = np.meshgrid(bins[:-1], range(len(algorithms)))
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros_like(xpos)
    
    dx = dy = 0.5
    dz = frequency_data.flatten()

    colors = plt.cm.plasma(np.linspace(0, 1, len(algorithms)))
    bar_colors = np.repeat(colors, len(bins) - 1, axis=0)
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color = bar_colors, alpha = 0.85, edgecolor='k')

    ax.set_xlabel("Scores", fontsize = 12)
    ax.set_ylabel("Algorithms", fontsize = 12)
    ax.set_yticks(range(len(algorithms)))
    ax.set_yticklabels(algorithms, fontsize = 10)
    ax.set_zlabel("Normalized Frequency", fontsize = 12)
    plt.title("3D Histogram of Normalized Score Distribution by Algorithm", fontsize = 14)
    plt.tight_layout()
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