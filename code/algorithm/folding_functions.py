import random as rd

def two_strings_fold(protein_sequence: str) -> list[int]:
    sequence_list: list[int] = []
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

def random_fold(protein_sequence: str) -> list[int]:
    sequence_list: list[int] = []
    protein_length = len(protein_sequence)
    random_choice = [-2, -1, 1, 2]
    previous_direction = 0

    for _ in range(protein_length - 1):
        direction = rd.choice(random_choice)
        while direction == previous_direction:
            direction = rd.choice(random_choice)

        sequence_list.append(direction)
        previous_direction = -1 * direction

    sequence_list.append(0)
    return sequence_list