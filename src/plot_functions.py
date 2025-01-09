import matplotlib.pyplot as plt
import numpy as np
from protein import Protein

def visualize(protein: Protein) -> None:
    Plot(protein)
    
class Plot():
    def __init__(self, protein: Protein) -> None:
        self.protein_sequence = protein.protein_sequence
        self.structure = protein.structure
        self.plot_structure()
        
    def plot_structure(self) -> None:
        plt.figure(figsize = (10, 10))

        color_map = {'H' : 'red', 'P' : 'blue', 'C' : 'green'}

        x_coords = [p[0] for p in self.structure]
        y_coords = [p[1] for p in self.structure]

        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)

        # Set limits and distance between amino acids
        plt.xlim(x_min + 5, x_max + 5)
        plt.ylim(y_min + 5, y_max + 5)
        plt.xticks(np.arange(x_min - 5, x_max + 5, 1))
        plt.yticks(np.arange(y_min - 5, y_max + 5, 1))

        # Plot amino acids as dots
        for i, (x, y) in enumerate(self.structure):
            amino_type = self.protein_sequence[i]
            plt.scatter(x, y, color = color_map[amino_type], s = 100, zorder = 3)

        # Plot sequential connections
        self._plot_sequential_connections()

        # Plot polar bonds
        self._plot_polar_connections()
            
        # Dummy plot for legend
        for amino_type, colour in color_map.items():
           plt.scatter([], [], color=colour, label = amino_type)

        plt.title('2D Protein Plot')
        plt.legend(loc = 'best')
        plt.axis('off')
        plt.show()

    def _plot_sequential_connections(self) -> None:
        x_prev, y_prev = 0, 0 
        for i, (x, y) in enumerate(self.structure):
            if i == 0:
                x_prev, y_prev = x, y
            else:
                plt.plot([x_prev, x], [y_prev, y], color = 'black', linestyle ='-', linewidth = 2, zorder= 2)
                x_prev, y_prev = x, y
    
    def _plot_polar_connections(self) -> None:
         for self.x_polar, self.y_polar in self.structure:
            if (self.structure[self.x_polar, self.y_polar][0] != 'P'):
                if (self.x_polar - 1, self.y_polar) in self.structure:
                    if self._check_neighbour(-1):
                        plt.plot([self.x_polar, self.x_polar - 1], [self.y_polar, self.y_polar], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar + 1, self.y_polar) in self.structure:        
                    if self._check_neighbour(1):
                        plt.plot([self.x_polar, self.x_polar + 1], [self.y_polar, self.y_polar], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar, self.y_polar - 1) in self.structure:        
                    if self._check_neighbour(-2):
                        plt.plot([self.x_polar, self.x_polar], [self.y_polar, self.y_polar - 1], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar, self.y_polar + 1) in self.structure:        
                    if self._check_neighbour(2):
                        plt.plot([self.x_polar, self.x_polar], [self.y_polar, self.y_polar + 1], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
    
    def _check_neighbour(self, direction: int) -> bool:
        if direction == -1 and not self._check_sequential(self.x_polar - 1, self.y_polar):
            dx, dy = -1, 0
        elif direction == 1 and not self._check_sequential(self.x_polar + 1, self.y_polar):
            dx, dy = 1, 0
        elif direction == -2 and not self._check_sequential(self.x_polar, self.y_polar - 1):
            dx, dy = 0, -1
        elif direction == 2 and not self._check_sequential(self.x_polar, self.y_polar + 1):
            dx, dy = 0, 1
        else:
            return False
        
        amino_1 = self.structure[self.x_polar, self.y_polar][0]
        amino_2 = self.structure[self.x_polar + dx, self.y_polar + dy][0]
        return self._check_pair(amino_1, amino_2)

    def _check_sequential(self, x: int, y: int) -> bool:
        i = self.structure[(self.x_polar, self.y_polar)][1]
        return self.structure[(x, y)][1] in [(i - 1), (i + 1)]

    def _check_pair(self, amino_1: str, amino_2: str) -> bool:
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return True
        elif amino_1 == 'C' and amino_2 == 'C':
            return True
        return False