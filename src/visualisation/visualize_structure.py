import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from src.classes import Protein
from typing import Optional

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
    
    protein_structure: dict[tuple[int, int, int], tuple[str, int]] = {}
    protein_sequence = protein.protein_sequence
    protein_rating = protein.get_rating()

    if protein.structure is not None:
        protein_structure = protein.structure.get_structure()

    x_coords = [p[0] for p in protein_structure]
    y_coords = [p[1] for p in protein_structure]
    z_coords = [p[2] for p in protein_structure]

    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)
    z_min, z_max = min(z_coords), max(z_coords)

    if dimension == 2:
        _plot_2d(x_min, y_min, x_max, y_max, protein_sequence, protein_structure, protein_rating, algorithm)

    elif dimension == 3: 
        _plot_3d(x_min, y_min, z_min, x_max, y_max, z_max, protein_sequence, protein_structure, protein_rating, algorithm)

    if save:
        if algorithm == "Brute Force":
            plt.savefig(f"{file_path}.png", dpi = 600)
        else:
            plt.savefig(f"{file_path}/3D {algorithm} Best Fold.png", dpi = 600)
    if show:
        plt.show()

    plt.close()

def _plot_2d(
    x_min: int,
    y_min: int,
    x_max: int,
    y_max: int,
    protein_sequence: str,
    protein_structure: dict[tuple[int, int, int], tuple[str, int]],
    protein_rating: int,
    algorithm: str,
    ) -> None:
    """Helper function for plotting the structure in 2D."""
    plt.figure(figsize = (12, 12))
    plt.title(f"2D protein plot\nAlgortihm: {algorithm} \nProtein: {protein_sequence}\nscore: {protein_rating}", fontsize = 12, fontweight = "bold")

    # Set limits and distance between amino acids
    plt.xlim(x_min - 5, x_max + 5)
    plt.ylim(y_min - 5, y_max + 5)
    plt.xticks(np.arange(x_min - 5, x_max + 5, 1))
    plt.yticks(np.arange(y_min - 5, y_max + 5, 1))
    plt.axis("off")

    # Plot amino acids as dots
    _plot_amino_acids(protein_structure)

    # Plot sequential connections
    _plot_sequential_connections(protein_structure)

    # Plot polar bonds
    _plot_polar_connections(protein_structure)

    # Dummy plot for legend
    _plot_legend()
    plt.legend(loc = "best")

def _plot_3d(
    x_min: int,
    y_min: int,
    z_min: int,
    x_max: int,
    y_max: int,
    z_max: int,
    protein_sequence: str,
    protein_structure: dict[tuple[int, int, int], tuple[str, int]],
    protein_rating: int,
    algorithm: str,
    ) -> None:
    """Helper function for plotting the structure in 3D."""
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_subplot(projection='3d')

    ax.set_title(f"3D protein plot\nAlgortihm: {algorithm} \nProtein: {protein_sequence}\nscore: {protein_rating}", fontsize = 12, fontweight = "bold")

    # Set limits and distance between amino acids
    ax.set_xlim(x_min - 2, x_max + 2)
    ax.set_ylim(y_min - 2, y_max + 2)
    ax.set_zlim(z_min - 2, z_max + 2)

    # Remove labels
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")
    ax.xaxis.set_tick_params(which = "both", labelbottom = False, bottom = False)
    ax.yaxis.set_tick_params(which = "both", labelleft = False, left = False)
    ax.zaxis.set_tick_params(which = "both", labelleft = False, left = False)

    ax.grid(True, color = "grey", linestyle = "--", alpha = 0.85)

    # Plot amino acids as dots
    _plot_amino_acids(protein_structure, ax)
    
    # Plot sequential connections
    _plot_sequential_connections(protein_structure, ax)

    # Plot polar bonds
    _plot_polar_connections(protein_structure, ax)

    # Plot legend
    _plot_legend(ax)
    ax.legend(loc = "best")

def _plot_amino_acids(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]],
    ax: Optional[Axes] = None
    ) -> None:
    """Helper function for plotting the amino acids."""
    color_map: dict[str, str] = {"H" : "red", "P" : "blue", "C" : "green"}

    for (x, y, z), (amino_type, _) in protein_structure.items():
        if ax:
            ax.scatter(x, y, z, color = color_map[amino_type], s = 100, zorder = 3)

        else:
            plt.scatter(x, y, color = color_map[amino_type], s = 100, zorder = 3)
    
def _plot_sequential_connections(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]],
    ax: Optional[Axes] = None
    ) -> None:
    """Helper function for plotting the sequential connections of the protein with a black line."""
    coordinates: list[tuple[int, int, int]] = list(protein_structure.keys())

    for (x_prev, y_prev, z_prev), (x, y, z) in zip(coordinates[:-1], coordinates[1:]):
        if ax:
            ax.plot([x_prev, x], [y_prev, y], [z_prev, z], color = "black", linestyle = "-", linewidth = 2, zorder = 2)

        else:
            plt.plot([x_prev, x], [y_prev, y], color = "black", linestyle = "-", linewidth = 2, zorder= 2)

def _plot_polar_connections(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]],
    ax: Optional[Axes] = None
    ) -> None:
    """
    Helper function for plotting the polar connections of the protein.
    Type of connections is highlighted by the color of the dashed line.
    lime green: H-H connection and C-H connection
    darkorange: C-C connection
    """
    neighbour_map: set[tuple[int, int, int]] = {(1, 0, 0) , (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)}

    for x_current, y_current, z_current in protein_structure:

        amino_1 = protein_structure[(x_current, y_current, z_current)][0]
        # Check for polar amino acid
        if amino_1 != "P":
            # Find neighbouring coordinates
            for (dx, dy, dz) in neighbour_map:
                x_next = x_current + dx
                y_next = y_current + dy
                z_next = z_current + dz

                # Check for valid and non sequential points
                if ((x_next, y_next, z_next) in protein_structure and not 
                    _check_sequential(protein_structure, x_current, y_current, x_next, y_next, z_current, z_next)):

                    amino_2 = protein_structure[(x_next, y_next, z_next)][0]
                    color = _check_connection_type(amino_1, amino_2)
                    if color is not None:
                        if ax:
                            ax.plot(
                                [x_current, x_next],
                                [y_current, y_next],
                                [z_current, z_next],
                                color = color,
                                linestyle = "--",
                                linewidth = 2,
                                zorder = 1)
                        else:
                            plt.plot(
                                [x_current, x_next],
                                [y_current, y_next],
                                color = color,
                                linestyle = "--",
                                linewidth = 2,
                                zorder = 1)

def _check_sequential(
    protein_structure: dict[tuple[int, int, int], tuple[str, int]], 
    x_old: int, 
    y_old: int, 
    x_new: int, 
    y_new: int,
    z_old: int,
    z_new: int
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

def _plot_legend(ax: Optional[Axes] = None) -> None:
    """
    Helper function for plotting the legend manually with dummy plots.
    It adds the labels for the colored amino acids.
    It also adds the type of connections between polar amino acids.
    """
    color_map: dict[str, str] = {"H" : "red", "P" : "blue", "C" : "green"}

    polar_lines: list[dict[str, str]] = [
    {"color" : "lime", "label" : "H-H Connection"},
    {"color" : "lime", "label" : "H-C Connection"},
    {"color" : "darkorange", "label" : "C-C Connection"}]

    # Amino colors
    for amino_type, color in color_map.items():
        if ax:
            ax.scatter([], [], [], color = color, label = amino_type, s = 100)
        else:
            plt.scatter([], [], color=color, label = amino_type)

    # Polar lines
    for line in polar_lines:
        if ax:
            ax.plot([], [], [], color=line["color"], linestyle = "--", label = line["label"], linewidth = 2)
        else:
            plt.plot([], [], color = line["color"], linestyle = "--", label = line["label"], linewidth = 2)