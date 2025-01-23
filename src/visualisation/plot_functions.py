import matplotlib.pyplot as plt
import numpy as np

from src.classes import Protein
from src.utils.constants import color_map
from typing import Optional

# ====================================================================================================
# Start Visualize protein structure
# ====================================================================================================
def visualize_protein(
    dimension: int,
    protein: Protein,
    algorithm: str,
    file_path: str,
    show: bool,
    save: bool
    ) -> None:
    """
    Main function for visualizing the protein structure.
    Amino acids are represented as colored dots,
    where 'H' is red, 'P' is blue, and 'C' is green.
    Sequential connections are black lines.
    Polar connections are connections between 'H' and 'H', 'H' and 'C' or 'C' and 'C' amino acids.
    These connections are highlighted by colored dashed lines.

    Parameters
    ----------
    dimension : int
        The dimension in which the folding takes place (`2` or `3`).

    protein : `Protein` object
        Protein object containing the sequence and folding structure.

    algorithm : str
        The name of the algorithm used (for example `Hill Climber`).
    
    file_path : str
        Path to save the plots.

    show : bool
        If `True` show the plots.

    save : bool
        If `True` save the plots.
    """

    if protein.protein_rating == 1:
        print("Error: did not get a valid protein for the plot.")
        return None

    if dimension == 2:
        _plot_2d(protein, algorithm, show, save, file_path)
    elif dimension == 3: 
        _plot_3d(protein, algorithm, show, save, file_path)

def _plot_2d(
    protein: Protein,
    algorithm: str, 
    show: bool, 
    save: bool, 
    file_path: str
    ) -> None:
    """Helper function for plotting the structure in 2D."""
    protein_sequence: str = protein.protein_sequence
    protein_structure: dict[tuple[int, int, int], tuple[str, int]] = {}

    if protein.structure is not None:
        protein_structure = protein.structure.get_structure()

    x_coords = [p[0] for p in protein_structure]
    y_coords = [p[1] for p in protein_structure]

    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)

    plt.figure(figsize = (12, 12))
    # Set limits and distance between amino acids
    plt.xlim(x_min - 5, x_max + 5)
    plt.ylim(y_min - 5, y_max + 5)
    plt.xticks(np.arange(x_min - 5, x_max + 5, 1))
    plt.yticks(np.arange(y_min - 5, y_max + 5, 1))

    # Plot amino acids as dots
    for i, (x, y, _) in enumerate(protein_structure):
        amino_type = protein_sequence[i]
        plt.scatter(x, y, color = color_map[amino_type], s = 100, zorder = 3)

    # Plot sequential connections
    _plot_sequential_connections(protein_structure)

    # Plot polar bonds
    _plot_polar_connections(protein_structure)
        
    # Dummy plot for legend
    _plot_legend()

    plt.title(f"2D protein plot\nAlgortihm: {algorithm} \nProtein: {protein_sequence}\nscore: {protein.get_rating()}", fontsize = 12, fontweight = "bold")
    plt.legend(loc = "best")
    plt.axis("off")

    if save:
        if algorithm == "Brute Force":
            plt.savefig(f"{file_path}.png", dpi = 600)
        else:
            plt.savefig(f"{file_path}/2D {algorithm} Best Fold.png", dpi = 600)
    if show:
        plt.show()

    plt.close()

def _plot_3d(
    protein: Protein,
    algorithm: str, 
    show: bool, 
    save: bool, 
    file_path: str
    ) -> None:
    """Helper function for plotting the structure in 3D."""
    protein_sequence: str = protein.protein_sequence
    protein_structure: dict[tuple[int, int, int], tuple[str, int]] = {}

    if protein.structure is not None:
        protein_structure = protein.structure.get_structure()

    x_coords = [p[0] for p in protein_structure]
    y_coords = [p[1] for p in protein_structure]
    z_coords = [p[2] for p in protein_structure]

    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)

    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_subplot(projection='3d')
    ax.set_title(f"3D protein plot\nAlgortihm: {algorithm} \nProtein: {protein_sequence}\nscore: {protein.get_rating()}", fontsize = 12, fontweight = "bold")

    # Set limits and distance between amino acids
    ax.set_xlim(x_min - 2, x_max + 2)
    ax.set_ylim(y_min - 2, y_max + 2)
    ax.set_zlim(z_min - 2, z_max + 2)

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    ax.xaxis.set_tick_params(which = "both", labelbottom = False, bottom = False)
    ax.yaxis.set_tick_params(which = "both", labelleft = False, left = False)
    ax.zaxis.set_tick_params(which = "both", labelleft = False, left = False)

    ax.grid(True, color = "grey", linestyle = "--", alpha = 0.85)

    for i, (x, y, z) in enumerate(protein_structure):
        amino_type = protein_sequence[i]
        ax.scatter(x, y, z, color = color_map[amino_type], s = 100, zorder = 3)
    
    # Plot sequential connections
    x_prev, y_prev, z_prev = 0, 0, 0
    for i, (x, y, z) in enumerate(protein_structure):
        if i == 0:
            x_prev, y_prev, z_prev = x, y, z
        else:
            ax.plot([x_prev, x], [y_prev, y], [z_prev, z], color = "black", linestyle = "-", linewidth = 2, zorder = 2)
            x_prev, y_prev, z_prev = x, y, z

    # Plot polar bonds
    direction_map_3d = {(1, 0, 0) , (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)}
    for x_current, y_current, z_current in protein_structure:

        amino_1 = protein_structure[(x_current, y_current, z_current)][0]
        # Check for polar amino acid
        if amino_1 != "P":
            # Find neighbouring coordinates
            for (dx, dy, dz) in direction_map_3d:
                x_next = x_current + dx
                y_next = y_current + dy
                z_next = z_current + dz

                # Check for valid and non sequential points
                if (x_next, y_next, z_next) in protein_structure and not _check_sequential(protein_structure, x_current, y_current, x_next, y_next, z_current, z_next):
                    amino_2 = protein_structure[(x_next, y_next, z_next)][0]

                    color = _check_connection_type(amino_1, amino_2)
                    if color is not None:
                        ax.plot([x_current, x_next], [y_current, y_next], [z_current, z_next], color = color, linestyle = "--", linewidth = 2, zorder = 1)
    # Plot legend
    for amino_type, color in color_map.items():
        ax.scatter([], [], [], color = color, label = amino_type, s = 100)
    # Dummy plots for polar connections
    polar_lines = [
        {"color" : "lime", "label" : "H-H Connection"},
        {"color" : "lime", "label" : "H-C Connection"},
        {"color" : "darkorange", "label" : "C-C Connection"}]
    for line in polar_lines:
        ax.plot([], [], [], color=line["color"], linestyle="--", label=line["label"], linewidth=2)

    ax.legend(loc = "best")

    if save:
        if algorithm == "Brute Force":
            plt.savefig(f"{file_path}.png", dpi = 600)
        else:
            plt.savefig(f"{file_path}/3D {algorithm} Best Fold.png", dpi = 600)
    if show:
        plt.show()

    plt.close()

def _plot_sequential_connections(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]]
    ) -> None:
    """Helper function for plotting the sequential connections of the protein with a black line."""
    x_prev, y_prev = 0, 0 
    for i, (x, y, _) in enumerate(protein_structure):
        if i == 0:
            x_prev, y_prev = x, y
        else:
            plt.plot([x_prev, x], [y_prev, y], color = "black", linestyle = "-", linewidth = 2, zorder= 2)
            x_prev, y_prev = x, y

def _plot_polar_connections(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]]
    ) -> None:
    """
    Helper function for plotting the polar connections of the protein.
    Type of connections is highlighted by the color of the dashed line.
    lime green: H-H connection and C-H connection
    darkorange: C-C connection
    """
    # Mapping for all neighbouring points
    direction_map_2d = {(1, 0) , (-1, 0), (0, 1), (0, -1)}
    for x_current, y_current, _ in protein_structure:

        amino_1 = protein_structure[(x_current, y_current, 0)][0]
        # Check for polar amino acid
        if amino_1 != "P":
            # Find neighbouring coordinates
            for (dx, dy) in direction_map_2d:
                x_next = x_current + dx
                y_next = y_current + dy

                # Check for valid and non sequential points
                if (x_next, y_next, _) in protein_structure and not _check_sequential(protein_structure, x_current, y_current, x_next, y_next):
                    amino_2 = protein_structure[(x_next, y_next, 0)][0]

                    color = _check_connection_type(amino_1, amino_2)
                    if color is not None:
                        plt.plot([x_current, x_next], [y_current, y_next], color = color, linestyle = "--", linewidth = 2, zorder = 1)

def _check_sequential(
    protein_structure: dict[tuple[int, int], tuple[str, int]], 
    x_old: int, 
    y_old: int, 
    x_new: int, 
    y_new: int,
    z_old: int = 0,
    z_new: int = 0
    ) -> bool:
    """
    Helper function for checking if two amino acids are in a sequential order. 
    Returns True if sequential, False otherwise.
    """
    i = protein_structure[(x_old, y_old, z_old)][1]
    return protein_structure[(x_new, y_new, z_new)][1] in [(i - 1), (i + 1)]

def _check_connection_type(
    amino_1: str,
    amino_2: str
    ) -> Optional[str]:
    """
    Helper function for checking the pair for possible polar connections.
    Returns the type of connection with color, None if no polar connection.
    """
    if (amino_1 == "H" and amino_2 in ["H", "C"]) or (amino_1 == "C" and amino_2 == "H"):
        return "lime"
    elif amino_1 == "C" and amino_2 == "C":
        return "darkorange"
    return None

def _plot_legend() -> None:
    """
    Helper function for plotting the legend manually with dummy plots.
    It adds the labels for the colored amino acids.
    It also adds the type of connections between polar amino acids.
    """
    # Label polar connections with dashed colored lines
    plt.plot([], [], color = "lime", linestyle = "--", label = "H-H Connection")
    plt.plot([], [], color = "lime", linestyle = "--", label = "H-C Connection")
    plt.plot([], [], color = "darkorange", linestyle = "--", label = "C-C Connection")

    for amino_type, colour in color_map.items():
        plt.scatter([], [], color=colour, label = amino_type)

# ========================================================================================================
# End Visualize protein structure
# ========================================================================================================

# ========================================================================================================
# Start Histogram protein distribution for specific algorithm
# ========================================================================================================
def histogram(
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
# ========================================================================================================
# End Histogram protein distribution for specific algorithm
# ========================================================================================================

# ========================================================================================================
# Start Hill Climber/Simulated Annealing score progression plot
# ========================================================================================================
def score_progression(
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
        Protein sequence (for example `HHHPPPHPCCP`).

    score_list : list[int]
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