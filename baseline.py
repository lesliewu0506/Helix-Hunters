import matplotlib.pyplot as plt

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

def rating(grid):
    score = 0
    # Define starting point to be in center
    x_point = (len(grid) // 2) - 1
    y_point = (len(grid) // 2) - 1

    for i in range(len(sequence)):
        current_amino = grid[y_point][x_point]
        print(current_amino)

        # Check neighbours
        if current_amino == 'H':
            if grid[y_point + 1][x_point] == 'H':
                score -= 1
            if grid[y_point - 1][x_point] == 'H':
                score -= 1
            if grid[y_point][x_point + 1] == 'H': 
                score -= 1   
            if grid[y_point][x_point - 1] == 'H':
                score -= 1
        
        elif current_amino == 'C':
            pass
        
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
    Protein1 = "HHPHHHPH"
    Protein2 = "HHPHHHPHPHHHPH"
    Protein3 = "HPHPPHHPHPPHPHHPPHPH"
    Protein4 = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    Protein5 = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    sequence = two_strings_fold(Protein1)
    grid = add_sequence_to_grid(Protein1, sequence)
    rating(grid)