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
        return self.grid

# Rating functions
def rating(grid, protein_sequence):
    return count_adjacent(grid) + correction_to_rating(protein_sequence)

def count_adjacent(grid):
    rows = len(grid)
    cols = len(grid[0])
    count = 0

    for row in range(rows):
        for col in range(cols):
            if col < cols - 1 and grid[row][col] == "C" and grid[row][col + 1] == "C":
                count -= 5
            elif col < cols - 1 and (grid[row][col] in ["H", "C"]) and (grid[row][col + 1] in ["H", "C"]):
                count -= 1
            if row < rows - 1 and grid[row][col] == "C" and grid[row + 1][col] == "C":
                count -= 5
            elif row < rows - 1 and (grid[row][col] in ["H", "C"]) and (grid[row + 1][col] in ["H", "C"]):
                count -= 1

    return count

def correction_to_rating(protein_sequence):
    correction = 0
    for i in range(1, len(protein_sequence)):
        current_amino = protein_sequence[i]
        previous_amino = protein_sequence[i - 1]

        if current_amino == "H":
            if previous_amino in ["H", "C"]:
                correction += 1

        elif current_amino == "C":
            if previous_amino == "H":
                correction += 1
            elif previous_amino == "C":
                correction += 5

    return correction

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
    # Protein1 = "HHPHHHPH"
    # Protein2 = "HHPHHHPHPHHHPH"
    # Protein3 = "HPHPPHHPHPPHPHHPPHPH"
    # Protein4 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    # Protein5 = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    # sequence = two_strings_fold(Protein1)
    # protein = "HHPHPPPPH"
    # sequence = [1,2,-1,-1,2,2,1,-2,0]
    # grid = add_sequence_to_grid(protein, sequence)
    # print(rating(grid, protein))
    protein_sequence = "HHPHPPPPH"
    protein_structure = [1,2,-1,-1,2,2,1,-2,0]
    grid = Grid(protein_sequence, protein_structure)
    grid.get_grid()