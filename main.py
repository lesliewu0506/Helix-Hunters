from src.brute_force.brute_force import generate_all_foldings, brute_force
from src.algorithms.randomise import Random as rd

protein_sequences = ["HHPHHHPHPHHHPH",
                     "HPHPPHHPHPPHPHHPPHPH",
                     "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
                     "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
                     "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
                     "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
                     "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
                     "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]

if __name__ == "__main__":
    # =========================== Brute Force ==========================
    # NOTE: This function is good only for the first protein sequence.
    # Running it for the others results in huge data files (>5GB). 
    # This is only for demonstration purposes.

    # # Generate all possible folding structures
    # generate_all_foldings(protein_sequences[0])
    # # Try all possible combinations
    # brute_force(protein_sequences[0], save = True)
    pass
    # =========================== Random ==========================
    random = rd(protein_sequences[0])
    random.run(show_plot = True, save_plot = True)

    print(f"Best score for random algorithm: {random.best_score}")