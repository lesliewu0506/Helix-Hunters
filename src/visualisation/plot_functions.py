import matplotlib.pyplot as plt
import numpy as np

from src.classes.protein import Protein
from typing import Optional

def visualize(protein: Protein, show: bool = True, save: bool = False, file_path: str = "output") -> None:
    """
    Main function for visualizing the protein structure.
    It uses the Protein class to retrieve a structure.
    It has optional arguments:
    - show: Shows the plot. Default is always True.
    - save: Saves the plot. Default is always False.
    - file_path: Saves the plot in specific directory.
    """
    if protein.protein_rating == 1:
        return None

    _plot_structure(protein, show, save, file_path)
    
def _plot_structure(protein: Protein, show: bool, save: bool, file_path: str) -> None:
    """
    Plots the entire protein structure with the amino acids,
    sequential connections and polar connections.
    Amino acids are colored dots; H: red; P: blue; C: green.
    Sequential connections are black lines.
    Polar connections are highlighted by colored dashed lines.
    """
    protein_sequence: str = protein.protein_sequence
    protein_structure: dict[tuple[int, int], tuple[str, int]] = protein.structure

    plt.figure(figsize = (10, 10))

    color_map = {'H' : 'red', 'P' : 'blue', 'C' : 'green'}

    x_coords = [p[0] for p in protein.structure]
    y_coords = [p[1] for p in protein.structure]

    x_min, x_max = min(x_coords), max(x_coords)
    y_min, y_max = min(y_coords), max(y_coords)

    # Set limits and distance between amino acids
    plt.xlim(x_min + 5, x_max + 5)
    plt.ylim(y_min + 5, y_max + 5)
    plt.xticks(np.arange(x_min - 5, x_max + 5, 1))
    plt.yticks(np.arange(y_min - 5, y_max + 5, 1))

    # Plot amino acids as dots
    for i, (x, y) in enumerate(protein.structure):
        amino_type = protein_sequence[i]
        plt.scatter(x, y, color = color_map[amino_type], s = 100, zorder = 3)

    # Plot sequential connections
    _plot_sequential_connections(protein_structure)

    # Plot polar bonds
    _plot_polar_connections(protein_structure)
        
    # Dummy plot for legend
    _plot_legend(color_map)

    plt.title(f"2D protein plot \nProtein:{protein_sequence}\nscore: {protein.get_rating()}")
    plt.legend(loc = 'best')
    plt.axis('off')
    if save:
        plt.savefig(f"{file_path}.png", dpi = 600)
    if show:
        plt.show()
    else:
        plt.close()

def _plot_sequential_connections(protein_structure: dict[tuple[int, int], tuple[str, int]]) -> None:
    """Plots the sequential connections of the protein with a black line."""
    x_prev, y_prev = 0, 0 
    for i, (x, y) in enumerate(protein_structure):
        if i == 0:
            x_prev, y_prev = x, y
        else:
            plt.plot([x_prev, x], [y_prev, y], color = 'black', linestyle ='-', linewidth = 2, zorder= 2)
            x_prev, y_prev = x, y

def _plot_polar_connections(protein_structure: dict[tuple[int, int], tuple[str, int]]) -> None:
    """
    Plots the polar connections of the protein.
    Type of connections is highlighted by the color of the dashed line.
    lime: H-H connection and C-H connection
    darkorange: C-C connection
    """
    # Mapping for all neighbouring points
    direction_map = {(1, 0) , (-1, 0), (0, 1), (0, -1)}

    for x_current, y_current in protein_structure:

        amino_1 = protein_structure[(x_current, y_current)][0]
        # Check for polar amino acid
        if amino_1 != 'P':
            # Find neighbouring coordinates
            for (dx, dy) in direction_map:
                x_next = x_current + dx
                y_next = y_current + dy

                # Check for valid and non sequential points
                if (x_next, y_next) in protein_structure and not _check_sequential(protein_structure, x_current, y_current, x_next, y_next):
                    amino_2 = protein_structure[(x_next, y_next)][0]

                    color = _check_connection_type(amino_1, amino_2)
                    _helper_plot_polar(x_current, y_current, x_next, y_next, color)

def _check_sequential(protein_structure: dict[tuple[int, int], tuple[str, int]], x_old: int, y_old: int, x_new: int, y_new: int) -> bool:
    """
    Checks if two amino acids are in a sequential order. 
    Returns True if sequential, False otherwise.
    """        
    i = protein_structure[(x_old, y_old)][1]
    return protein_structure[(x_new, y_new)][1] in [(i - 1), (i + 1)]

def _check_connection_type(amino_1: str, amino_2: str) -> Optional[str]:
    """Checks the pair for possible connections and returns the type of connection with color."""
    if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
        return "lime"
    elif amino_1 == 'C' and amino_2 == 'C':
        return "darkorange"
    return None

def _helper_plot_polar(x_old: int, y_old: int, x_new: int, y_new: int, colour: Optional[str]) -> None:
    """Helper functions that plots the polar connections."""
    if colour is not None:
        plt.plot([x_old, x_new], [y_old, y_new], color = colour, linestyle = '--', linewidth = 2, zorder = 1)

def _plot_legend(color_map: dict[str, str]) -> None:
    """
    Plots the legend manually with dummy plots.
    It adds the labels for the colored amino acids.
    It also adds the type of connections between polar amino acids.
    """
    # Label polar connections with dashed colored lines
    plt.plot([], [], color = "lime", linestyle = '--', label = "H-H Connection")
    plt.plot([], [], color = "lime", linestyle = '--', label = "H-C Connection")
    plt.plot([], [], color = "darkorange", linestyle = '--', label = "C-C Connection")

    for amino_type, colour in color_map.items():
        plt.scatter([], [], color=colour, label = amino_type)

def histogram(protein_sequence: str, score: list[int], iterations: int, show: bool = False, save: bool = False, file_path: str = "output", algorithm: str = None) -> None:
    """Creates a stylish histogram with gradient color and improved aesthetics."""
    plt.figure(figsize=(12, 7))

    # Center bars on each x-tick
    bins = np.arange(min(score) - 0.5, max(score) + 1.5, 1)
    x_ticks = np.arange(min(score), max(score) + 1, 1)

    n, bins, patches = plt.hist(score, bins = bins, edgecolor = 'black', alpha = 0.85, linewidth = 1.5)
    
    # Set different color for bars
    cmap = plt.get_cmap('plasma')
    for patch, value in zip(patches, n):
        patch.set_facecolor(cmap(value / max(n)))
    
    plt.title(f'Protein Score Distribution\nProtein sequence:{protein_sequence}\nAlgorithm:{algorithm}, {iterations} iterations', fontsize=12, fontweight='bold', pad=20)
    plt.xlabel('Protein Score', fontsize=14, labelpad=15)
    plt.ylabel('Frequency', fontsize=14, labelpad=15)
    
    plt.xticks(x_ticks, fontsize=12)
    plt.yticks(fontsize=12)
    
    # Add y-grid for clarity
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    
    plt.tight_layout()

    if save:
        string = f"Protein_score_distribution_{iterations}_iterations"
        plt.savefig(f"{file_path}/{string}.png", dpi = 600)

    if show:
        plt.show()
    else:
        plt.close()