# import matplotlib.pyplot as plt

def initialize_grid(protein_sequence):
    # Create 2D array with length 2 times sequence length plus 2 buffer
    grid_length = 2 * len(protein_sequence) + 2
    array = []

    for _ in range(grid_length):
        column = []
        for _ in range(grid_length):
            column.append(0)
        array.append(column)

    return array

def add_sequence_to_grid(protein_sequence, sequence):
    # Initialize grid
    grid = initialize_grid(protein_sequence)

    # Define starting point to be in center
    x_point = (len(grid) // 2) - 1
    y_point = (len(grid) // 2) - 1

    # Add points to grid
    for i in range(len(sequence)):
        grid[y_point][x_point] = protein_sequence[i]
        x_point, y_point = update_position(sequence, i, x_point, y_point) 
    
    # Print grid
    for i in range(len(grid)):
        print(grid[i])
    return grid

def update_position(sequence, i, x_coord, y_coord):
    if sequence[i] == 1:
        x_coord += 1
    elif sequence[i] == -1:
        x_coord -= 1
    elif sequence[i] == 2:
        y_coord -= 1
    elif sequence[i] == -2:
        y_coord += 1
    return x_coord, y_coord

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
    Protein1 = "HHPHHHPH"
    Protein2 = "HHPHHHPHPHHHPH"
    Protein3 = "HPHPPHHPHPPHPHHPPHPH"
    Protein4 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    Protein5 = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    sequence = two_strings_fold(Protein1)
    protein = "HHPHPPPPH"
    sequence = [1,2,-1,-1,2,2,1,-2,0]
    grid = add_sequence_to_grid(protein, sequence)
    print(rating(grid, protein))