import csv
import matplotlib.pyplot as plt
import numpy as np

class Protein():
    def __init__(self, protein_sequence, function):
        self.protein_sequence = protein_sequence
        self._add_folding_structure(function)

    def _add_folding_structure(self, function):
        self.amino_direction = function(self.protein_sequence)

        self.structure = Grid(self.protein_sequence, self.amino_direction).get_structure()
        self.protein_rating = Rating(self.protein_sequence, self.structure)

    def output_csv(self):
        with open('output.csv', 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['amino', 'fold'])
            for amino, direction in zip(self.protein_sequence, self.amino_direction):
                writer.writerow([amino, direction])
            writer.writerow(['score', self.protein_rating.score])

class Grid():
    def __init__(self, protein_sequence, amino_direction):
        self.structure = {}
        self._create_structure(protein_sequence, amino_direction)

    def _create_structure(self, protein_sequence, amino_direction):
        self.x_coord = 0 
        self.y_coord = 0

        for i in range(len(protein_sequence)):
            self.structure[(self.x_coord, self.y_coord)] = (protein_sequence[i], i)
            self._update_position(amino_direction[i])

    def _update_position(self, direction):
        if direction == 1:
            self.x_coord += 1
        elif direction == -1:
            self.x_coord -= 1
        elif direction == 2:
            self.y_coord += 1
        elif direction == -2:
            self.y_coord -= 1

    def get_structure(self):
        return self.structure

class Rating():
    def __init__(self, protein_sequence, structure: dict):
        self.protein_sequence = protein_sequence
        self.structure = structure

        self.score = 0
        self._count_adjacent()
        # Prevent double counting
        self.score = self.score // 2

    def _count_adjacent(self):
        direction_map = {1 : (1, 0), -1 : (-1, 0), 2 : (0, 1), -2 : (0, -1)}

        for self.x_current, self.y_current in self.structure:
            for direction, (dx, dy) in direction_map.items():
                x_next = self.x_current + dx
                y_next = self.y_current + dy

                if (x_next, y_next) in self.structure:
                    self.score += self._check_neighbour(direction)

    def _check_neighbour(self, direction):
        if direction == -1 and not self._check_sequential(self.x_current - 1, self.y_current):
            dx, dy = -1, 0
        elif direction == 1 and not self._check_sequential(self.x_current + 1, self.y_current):
            dx, dy = 1, 0
        elif direction == -2 and not self._check_sequential(self.x_current, self.y_current - 1):
            dx, dy = 0, -1
        elif direction == 2 and not self._check_sequential(self.x_current, self.y_current + 1):
            dx, dy = 0, 1
        else:
            dx, dy = 0, 0
            return 0
        
        amino_1 = self.structure[self.x_current, self.y_current][0]
        amino_2 = self.structure[self.x_current + dx, self.y_current + dy][0]
        return self._check_pair(amino_1, amino_2)

    def _check_sequential(self, x, y):
        i = self.structure[(self.x_current, self.y_current)][1]
        return self.structure[(x, y)][1] in [(i - 1), (i + 1)]

    def _check_pair(self, amino_1, amino_2):
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return -1
        elif amino_1 == 'C' and amino_2 == 'C':
            return -5
        else:
            return 0
    
    def get_rating(self):
        print(self.score)
        return self.score

class Plot():
    def __init__(self, protein: Protein):
        self.protein_sequence = protein.protein_sequence
        self.structure = protein.structure
    
    def plot_structure(self):
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

    def _plot_sequential_connections(self):
        x_prev, y_prev = 0, 0 
        for i, (x, y) in enumerate(self.structure):
            if i == 0:
                x_prev, y_prev = x, y
            else:
                plt.plot([x_prev, x], [y_prev, y], color = 'black', linestyle ='-', linewidth = 2, zorder= 2)
                x_prev, y_prev = x, y
    
    def _plot_polar_connections(self):
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
    
    def _check_neighbour(self, direction):
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

    def _check_sequential(self, x, y):
        i = self.structure[(self.x_polar, self.y_polar)][1]
        return self.structure[(x, y)][1] in [(i - 1), (i + 1)]

    def _check_pair(self, amino_1, amino_2):
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return True
        elif amino_1 == 'C' and amino_2 == 'C':
            return True

def two_strings_fold(protein_sequence):
    sequence_list = []
    protein_length = len(protein_sequence)
    for i in range(protein_length):
        if i < (protein_length // 2 - 1):
            sequence_list.append(1)
        elif i == (protein_length // 2 - 1):
            sequence_list.append(2)
        elif i == (protein_length - 1):
            sequence_list.append(0)
        else:
            sequence_list.append(-1)
    return sequence_list

# def main(protein_sequence, function):
    # protein = Protein(protein_sequence)
    # protein.add_folding_structure(function)
    # protein.protein_grid.get_grid()
    # protein.protein_rating.get_rating()
    # protein.output_csv()

if __name__ == "__main__":
    protein_sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
    # main(protein_sequence, two_strings_fold)
    protein = Protein(protein_sequence, two_strings_fold)
    protein.protein_rating.get_rating()
    plot = Plot(protein)
    plot.plot_structure()