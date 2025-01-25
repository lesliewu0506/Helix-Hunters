"""
This module contains general-purpose constants that can be used
across the project. These constant are designed to simplify repetitive tasks, improve 
code readability, and promote reusability.

Constants in this file are not specific to any single module but can be 
used by various parts of the project wherever needed.
"""
# ===============================================================
# Greedy Constants
# ===============================================================

DIRECTIONS_2D: list[int] = [0, 1, 2]
DIRECTIONS_3D: list[int] = [0, 1, 2, 3, 4] 

# ===============================================================
# Hill Climber Constants
# ===============================================================

MAX_UNCHANGED_ITERATIONS: int = 600
DIRECTION_CHOICES_2D: list[int] = [-2, -1, 1, 2]
DIRECTION_CHOICES_3D: list[int] = [-3, -2, -1, 1, 2, 3]

# ===============================================================
# Simulated Annealing Constants
# ===============================================================

TEMPERATURE: int = 4
TEMPERATURE_DECAY: int = 0.999

# ===============================================================
# Utils Constant
# ===============================================================

ITERATIVE_ALGORITHM_FACTOR: int = 10

# ===============================================================
# Protein Sequences
# ===============================================================

PROTEIN_SEQUENCES: list[str] = [
    "HHPHHHPH",
    "HHPHHHPHPHHHPH",
    "HPHPPHHPHPPHPHHPPHPH",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]

# ===============================================================
# Global Mappings
# ===============================================================

PROTEIN_SEQUENCE_MAP: dict[str, str] = {
    "HHPHHHPH" : "0",
    "HHPHHHPHPHHHPH" : "1",
    "HPHPPHHPHPPHPHHPPHPH" : "2",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

ALGORITHM_FOLDER_MAP: dict[str, str] = {
    "Random" : "random",
    "Greedy" : "greedy",
    "Hill Climber" : "hill",
    "Simulated Annealing" : "annealing",
    "Plant Propagation" : "propagation",
    "Genetic Algorithm" : "genetic"}

DIRECTION_MAP_2D: dict[int, list[int]] = {
    1 : [2, 1, -2],
    -1 : [-2, -1, 2],
    2 : [-1, 2, 1],
    -2 : [1, -2, -1]}

DIRECTION_MAP_3D: dict[int, list[int]] = {
    1 : [2, 1, -2, 3, -3],
    -1 : [-2, -1, 2, -3, 3],
    2 : [-1, 2, 1, 3, -3],
    -2 : [1, -2, -1, -3, 3],
    3 : [1, 3, -1, 2, -2],
    -3 : [-1, -3, 1, -2, 2]}

# ===============================================================
# Algorithms
# ===============================================================

ALGORTIHMS: list[str] = ["Random", "Greedy", "Hill Climber", "Simulated Annealing"]