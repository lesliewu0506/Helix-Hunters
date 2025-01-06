def initialize_grid(protein_sequence):
    # Create 2D array with length 2 times sequence length
    grid_length = 2 * len(protein_sequence)
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
        if sequence[i] == 1:
            x_point += 1
        elif sequence[i] == -1:
            x_point -= 1
        elif sequence[i] == 2:
            y_point -= 1
        elif sequence[i] == -2:
            y_point += 1
    
    # Print grid
    for i in range(len(grid)):
        print(grid[i])
    return

def rating():

    return

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
    Protein = "HHPHHHPH"
    sequence = two_strings_fold(Protein)
    add_sequence_to_grid(Protein, sequence)