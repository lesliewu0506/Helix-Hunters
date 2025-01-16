from src.brute_force.brute_force import generate_all_foldings, brute_force
from src.algorithms.randomise import Random as rd
from src.algorithms.greedy import Greedy
from src.algorithms.hill_climber import HillClimber

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

    # =========================== Random ===============================
    # random = rd(protein_sequences[0])
    # random.run(show_plot=False, save_data=True, repeats= 10, iterations = 10000)
    # for protein_sequence in protein_sequences:
    #     random = rd(protein_sequence)
    #     random.run(show_plot = False, save_plot = True)

    #     print(f"Best score for random algorithm for protein {protein_sequence}: {random.best_score}")

    # ====================== Random Greedy ==============================
    # greedy = Greedy(protein_sequences[0])
    # greedy.run(save_data=True, repeats = 10, iterations= 10000)

    # ====================== Hill Climber ===============================
    hillclimber = HillClimber(protein_sequences[6])
    hillclimber.run(iterations=3000)