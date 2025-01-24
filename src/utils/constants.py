"""
This module contains general-purpose constants that can be used
across the project. These constant are designed to simplify repetitive tasks, improve 
code readability, and promote reusability.

Constants in this file are not specific to any single module but can be 
used by various parts of the project wherever needed.
"""
# ===============================================================
# Protein Sequences
# ===============================================================

protein_sequences: list[str] = [
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

protein_sequence_map: dict[str, str] = {
    "HHPHHHPH" : "0",
    "HHPHHHPHPHHHPH" : "1",
    "HPHPPHHPHPPHPHHPPHPH" : "2",
    "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP" : "3",
    "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH" : "4",
    "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP" : "5",
    "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC" : "6",
    "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH" : "7",
    "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH" : "8"}

algorithm_folder_map: dict[str, str] = {
    "Random" : "random",
    "Greedy" : "greedy",
    "Hill Climber" : "hill",
    "Simulated Annealing" : "annealing",
    "Plant Propagation" : "propagation",
    "Genetic Algorithm" : "genetic"}

direction_map_2d: dict[int, list[int]] = {
    1 : [2, 1, -2],
    -1 : [-2, -1, 2],
    2 : [-1, 2, 1],
    -2 : [1, -2, -1]}

direction_map_3d: dict[int, list[int]] = {
    1 : [2, 1, -2, 3, -3],
    -1 : [-2, -1, 2, -3, 3],
    2 : [-1, 2, 1, 3, -3],
    -2 : [1, -2, -1, -3, 3],
    3 : [1, 3, -1, 2, -2],
    -3 : [-1, -3, 1, -2, 2]}

# ===============================================================
# Algorithms
# ===============================================================

algorithms: list[str] = ["Random", "Greedy", "Hill Climber", "Simulated Annealing"]

