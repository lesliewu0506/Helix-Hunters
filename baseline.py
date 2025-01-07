import csv

class Protein():
    def __init__(self, protein_sequence):
        self.protein_sequence = protein_sequence

    def add_folding_structure(self, function):
        self.protein_structure = function(self.protein_sequence)
        self.function_name = function.__name__

        self.protein_grid = Grid(self.protein_sequence, self.protein_structure)
        self.protein_rating = Rating(self.protein_sequence, self.protein_structure, self.protein_grid.grid)

    def output_csv(self):
        with open('output.csv', 'w', newline = '') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['amino', 'fold'])
            for amino, direction in zip(self.protein_sequence, self.protein_structure):
                writer.writerow([amino, direction])
            writer.writerow(['score', self.protein_rating.score])

class Grid():
    def __init__(self, protein_sequence, protein_structure):

        self.protein_sequence = protein_sequence
        self.protein_structure = protein_structure

        # Create 2D array with length 2 times sequence length plus 2 buffer
        self.grid_size = 2 * len(self.protein_sequence) + 2

        self._initialize_grid()
        self._add_structure_to_grid()

    def _initialize_grid(self):
        self.grid = []
        for _ in range(self.grid_size):
            column = []
            for _ in range(self.grid_size):
                column.append(0)
            self.grid.append(column)

    def _add_structure_to_grid(self):
        # Define starting point to be in center
        self.x_coord = self.grid_size // 2
        self.y_coord = self.grid_size // 2

        # Add points to grid
        for i in range(len(self.protein_sequence)):
            self.grid[self.y_coord][self.x_coord] = self.protein_sequence[i]
            self._update_position(self.protein_structure[i]) 

    def _update_position(self, direction):
        if direction == 1:
            self.x_coord += 1
        elif direction == -1:
            self.x_coord -= 1
        elif direction == 2:
            self.y_coord -= 1
        elif direction == -2:
            self.y_coord += 1
    
    def get_grid(self):
        for i in range(self.grid_size):
            print(self.grid[i])

class Rating():
    def __init__(self, protein_sequence, protein_structure, grid):
        self.protein_sequence = protein_sequence
        self.protein_structure = protein_structure
        self.grid = grid

        self.score = 0

        self._count_adjacent()
        self._correction_to_rating()

    def _count_adjacent(self):
        rows = len(self.grid)
        cols = len(self.grid[0])

        for row in range(rows):
            for col in range(cols):
                if col < cols - 1 and self.grid[row][col] == "C" and self.grid[row][col + 1] == "C":
                    self.score -= 5
                elif col < cols - 1 and (self.grid[row][col] in ["H", "C"]) and (self.grid[row][col + 1] in ["H", "C"]):
                    self.score -= 1
                if row < rows - 1 and self.grid[row][col] == "C" and self.grid[row + 1][col] == "C":
                    self.score -= 5
                elif row < rows - 1 and (self.grid[row][col] in ["H", "C"]) and (self.grid[row + 1][col] in ["H", "C"]):
                    self.score -= 1

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

if __name__ == "__main__":
    protein_sequence = "HHPHPPPPH"
    # protein_structure = [1,2,-1,-1,2,2,1,-2,0]
    # grid = Grid(protein_sequence, protein_structure)
    # rating = Rating(protein_sequence, protein_structure, grid.grid)
    # print(rating.score)
    protein = Protein(protein_sequence)
    protein.add_folding_structure(two_strings_fold)
    protein.protein_grid.get_grid()
    protein.output_csv()
