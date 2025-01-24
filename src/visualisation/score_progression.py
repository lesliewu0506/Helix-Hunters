import numpy as np
import matplotlib.pyplot as plt

def visualize_score_progression(
    dimension: int,
    protein_sequence: str, 
    score_list: list[int],
    algorithm: str,
    file_path: str,
    show: bool,
    save: bool
    ) -> None:
    """
    Creates a plot for the evolution of the protein score with the Hill Climber/Simulated Annealing algorithm.
    
    Parameters
    ----------
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    protein_sequence : str
        Protein sequence (for example `HHPHHHPH`).

    score_list : `list[int]`
        List of the scores during a run.

    algorithm : str
        The name of the algorithm used (either `Hill Climber` or `Simulated Annealing`).

    file_path : str
        Path to save the plots.
    
    show : bool
        If `True` show the plot.

    save : bool
        If `True` save the plot.
    """
    plt.figure(figsize = (12, 7))
    plt.plot(np.arange(1, len(score_list) + 1), score_list, label = "Score Progression", linewidth = 2)
    
    plt.title(f"{dimension}D {algorithm} Progression\nProtein: {protein_sequence}", fontsize = 14, fontweight = "bold")
    plt.xlabel("Iterations", fontsize = 12)
    plt.ylabel("Score", fontsize = 12)

    # Set axis attributes
    ymin = min(score_list)
    ymax = max(score_list)

    plt.xlim(0, len(score_list))
    plt.ylim(ymin - 1, ymax + 1)
    plt.yticks(np.arange(ymin, ymax + 1, 1))

    plt.grid(axis = "y", linestyle = "--", alpha = 0.85)
    plt.legend(loc = "best")

    plt.tight_layout()

    if save:
        plt.savefig(f"{file_path}/{dimension}D {algorithm} Score Progression.png", dpi = 600)

    if show:
        plt.show()

    plt.close()