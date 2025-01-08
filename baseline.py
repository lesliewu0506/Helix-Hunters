import csv
import matplotlib.pyplot as plt
import numpy as np

class Protein():
    def __init__(self, protein_sequence, function):
        self.protein_sequence = protein_sequence
        self._add_folding_structure(function)

    def _add_folding_structure(self, function):
        self.amino_direction = function(self.protein_sequence)
        grid = Grid(self.protein_sequence, self.amino_direction)
        self.structure = grid.get_structure()
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
        self.structure = self._create_structure(protein_sequence, amino_direction)

    def _create_structure(self, protein_sequence, amino_direction):
        self.coordinates = {}
        self.x_coord = 0 
        self.y_coord = 0

        for i in range(len(protein_sequence)):
            self.coordinates[(self.x_coord, self.y_coord)] = (protein_sequence[i], i)
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
        print(self.structure)
        return self.structure

class Rating():
    def __init__(self, protein_sequence, structure: dict):
        self.protein_sequence = protein_sequence
        self.structure = structure

        self.score = 0

        self._count_adjacent()
        self._correction_to_rating()

    def _count_adjacent(self):
        for self.x_current, self.y_current in self.structure:
            if (self.x_current - 1, self.y_current) in self.structure:
                self.score += self._check_neighbour(-1)

            if (self.x_current + 1, self.y_current) in self.structure:        
                self.score += self._check_neighbour(1)

            if (self.x_current, self.y_current - 1) in self.structure:        
                self.score +=  self._check_neighbour(-2)

            if (self.x_current, self.y_current + 1) in self.structure:        
                self.score +=  self._check_neighbour(2)

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
        
        amino_1 = self.structure[self.x_current, self.y_current]
        amino_2 = self.structure[self.x_current + dx, self.y_current + dy]
        return self._check_pair(amino_1, amino_2)

    def _check_sequential(self, x, y):
        i = self.structure[(self.x_current, self.y_current)][1]
        return self.structure[(x, y)][1] in [(i - 1), (i + 1)]

    def _check_pair(self, amino_1, amino_2):
        if (amino_1 == 'H' and amino_2 in ['H', 'C']) or (amino_1 == 'C' and amino_2 == 'H'):
            return -1
        elif amino_1 == 'C' and amino_2 == 'C':
            return -5
    
    def _correction_to_rating(self):
        for i in range(1, len(self.protein_sequence)):
            current_amino = self.protein_sequence[i]
            previous_amino = self.protein_sequence[i - 1]

            if current_amino == "H":
                if previous_amino in ["H", "C"]:
                    self.score += 1

            elif current_amino == "C":
                if previous_amino == "H":
                    self.score += 1
                elif previous_amino == "C":
                    self.score += 5

    def get_rating(self):
        print(self.score)

class Plot():
    def __init__(self, protein: Protein):
        self.protein_sequence = protein.protein_sequence
        self.protein_structure = protein.structure

        self._assign_coordinates()

    def _assign_coordinates(self):
        self.coordinates = {}
        self.coordinates_list = []
        self.x_coord = 0 
        self.y_coord = 0

        for i in range(len(self.protein_sequence)):
            self.coordinates[(self.x_coord, self.y_coord)] = self.protein_sequence[i]
            self.coordinates_list.append((self.x_coord, self.y_coord))
            self._update_position_plot(self.protein_structure[i])

    def _update_position_plot(self, direction):
        if direction == 1:
            self.x_coord += 1
        elif direction == -1:
            self.x_coord -= 1
        elif direction == 2:
            self.y_coord += 1
        elif direction == -2:
            self.y_coord -= 1
    
    def plot_structure(self):
        plt.figure(figsize = (10, 10))

        color_map = {'H' : 'red', 'P' : 'blue', 'C' : 'green'}

        x_coords = [x[0] for x in self.coordinates]
        y_coords = [y[1] for y in self.coordinates]
        
        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)

        # Set limits and distance between amino acids
        plt.xlim(x_min + 5, x_max + 5)
        plt.ylim(y_min + 5, y_max + 5)
        plt.xticks(np.arange(x_min - 5, x_max + 5, 1))
        plt.yticks(np.arange(y_min - 5, y_max + 5, 1))

        # Plot amino acids as dots
        for i, (x, y) in enumerate(self.coordinates):
            amino_type = self.protein_sequence[i]
            plt.scatter(x, y, color = color_map[amino_type], s = 100, zorder = 3)

        # Plot sequential connections
        for i in range(len(self.coordinates) - 1):
            x_prev, y_prev = self.coordinates_list[i]
            x_next, y_next = self.coordinates_list[i + 1]
            plt.plot([x_prev, x_next], [y_prev, y_next], color = 'black', linestyle ='-', linewidth = 2, zorder= 2)

        # Plot polar bonds
        for self.x_polar, self.y_polar in self.coordinates:
            if (self.coordinates[self.x_polar, self.y_polar] != 'P'):
                if (self.x_polar - 1, self.y_polar) in self.coordinates:
                    if self._check_neighbour(-1):
                        plt.plot([self.x_polar, self.x_polar - 1], [self.y_polar, self.y_polar], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar + 1, self.y_polar) in self.coordinates:        
                    if self._check_neighbour(1):
                        plt.plot([self.x_polar, self.x_polar + 1], [self.y_polar, self.y_polar], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar, self.y_polar - 1) in self.coordinates:        
                    if self._check_neighbour(-2):
                        plt.plot([self.x_polar, self.x_polar], [self.y_polar, self.y_polar - 1], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
                if (self.x_polar, self.y_polar + 1) in self.coordinates:        
                    if self._check_neighbour(2):
                        plt.plot([self.x_polar, self.x_polar], [self.y_polar, self.y_polar + 1], color = 'grey', linestyle = '--', linewidth = 2, zorder = 1)
            
        # Dummy plot for legend
        for amino_type, colour in color_map.items():
           plt.scatter([], [], color=colour, label = amino_type)

        plt.title('2D Protein Plot')
        plt.legend(loc = 'best')
        plt.axis('off')
        plt.show()

    def _check_neighbour(self, direction):
        if direction == -1:
            amino_1 = self.coordinates[self.x_polar, self.y_polar]
            amino_2 = self.coordinates[self.x_polar - 1, self.y_polar]
            return self._check_pair(amino_1, amino_2)
        elif direction == 1:
            amino_1 = self.coordinates[self.x_polar, self.y_polar]
            amino_2 = self.coordinates[self.x_polar + 1, self.y_polar]
            return self._check_pair(amino_1, amino_2)
        if direction == -2:
            amino_1 = self.coordinates[self.x_polar, self.y_polar]
            amino_2 = self.coordinates[self.x_polar, self.y_polar - 1]
            return self._check_pair(amino_1, amino_2)
        if direction == 2:
            amino_1 = self.coordinates[self.x_polar, self.y_polar]
            amino_2 = self.coordinates[self.x_polar, self.y_polar + 1]
            return self._check_pair(amino_1, amino_2)

    def _check_pair(self, amino_1, amino_2):
        if amino_2 in ['H', 'C']:
            return True

    def get_coordinates(self):
        print(self.coordinates)

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
    protein_sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    # main(protein_sequence, two_strings_fold)
    protein = Protein(protein_sequence, two_strings_fold)
    # plot = Plot(protein)
    # plot.get_coordinates()
    # plot.plot_structure()