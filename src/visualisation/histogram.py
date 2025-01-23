import numpy as np
import matplotlib.pyplot as plt

def visualize_histogram(
    dimension: int,
    protein_sequence: str,
    score_list: list[int],
    iterations: int,
    algorithm: str,
    file_path: str,
    show: bool,
    save: bool
    ) -> None:
    """
    Creates a stylish histogram with gradient color for the score distribution of an algorithm.
    
    Parameters
    ----------
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    protein_sequence : str
        Protein sequence (for example `HHHPPPHPCCP`).

    score_list : list[int]
        List of the scores during the run.

    iterations : int
        Amount of iterations for the run.

    algorithm : str
        The name of the algorithm used (for example `Hill Climber`).

    file_path : str
        Path to save the plots.

    show : bool
        If `True` show the plot.

    save : bool
        If `True` save the plot.
    """
    plt.figure(figsize = (12, 7))

    min_score = min(score_list)
    max_score = max(score_list)
    # Center bars on each x-tick
    bins = np.arange(min_score - 0.5, max_score + 1.5, 1)

    n, bins, patches = plt.hist(score_list, bins = bins, edgecolor = "black", alpha = 0.85, linewidth = 1.5)

    # Set different color for bars
    cmap = plt.get_cmap("plasma")
    for patch, value in zip(patches, n):
        patch.set_facecolor(cmap(value / max(n)))
    
    plt.title(f"{dimension}D Protein Score Distribution\nProtein sequence: {protein_sequence}\nAlgorithm: {algorithm}, {iterations} iterations", fontsize = 14, fontweight = "bold")
    plt.xlabel("Protein Score", fontsize = 12)
    plt.ylabel("Frequency", fontsize = 12)

    plt.xticks(np.arange(min_score, max_score + 1, 1))

    # Adjust y limit to round value
    multiple = iterations / 100
    y_max = np.ceil(max(n) / multiple) * multiple
    plt.ylim(top = y_max)

    plt.grid(axis = "y", linestyle = "--", alpha = 0.85)
    
    plt.tight_layout()

    if save:
        string = f"{dimension}D Protein Score Distribution {algorithm}"
        plt.savefig(f"{file_path}/{string}.png", dpi = 600)

    if show:
        plt.show()
    
    plt.close()