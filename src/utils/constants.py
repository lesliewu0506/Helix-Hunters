# ===============================================================
# Protein Sequences
# ===============================================================

protein_sequences: list[str] = [
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

color_map: dict[str, str] = {"H" : "red", "P" : "blue", "C" : "green"}

# ===============================================================
# Algorithms
# ===============================================================

algorithms: list[str] = ["Random", "Greedy", "Hill Climber", "Simulated Annealing"]