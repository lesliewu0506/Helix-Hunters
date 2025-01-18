from src.brute_force.Brute_Force import brute_force
from src.brute_force.generate_all_folding_structures import generate_all_foldings, _check_valid_folding
from src.algorithms.randomise import Random 
from src.algorithms.greedy import Greedy
from src.algorithms.hill_climber import HillClimber
from src.algorithms.simulated_annealing import SimulatedAnnealing
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

    # Generate all possible folding structures
    # generate_all_foldings(protein_sequences[0])
    # # Try all possible combinations
    # brute_force(protein_sequences[0], save = True)

    # =========================== Random ===============================
    # for protein_sequence in protein_sequences:
    #     random = Random(protein_sequence)
    #     random.run(save_plot = True, save_data = True, repeats = 1)

    #     print(f"Best score for Random algorithm for protein {protein_sequence}: {random.best_score}")

    # ====================== Random Greedy ==============================
    # for protein_sequence in protein_sequences:
    #     greedy = Greedy(protein_sequence)
    #     greedy.run(save_plot = True, save_data = True, repeats = 10)

    #     print(f"Best score for Greedy algorithm for protein {protein_sequence}: {greedy.best_score}")
    # ====================== Hill Climber ===============================
    # for protein_sequence in protein_sequences:
    #     hillclimber = HillClimber(protein_sequence)
    #     hillclimber.run(save_plot = True, save_data= True, repeats = 10)

    #     print(f"Best score for Hill Climber algorithm for protein {protein_sequence}: {hillclimber.best_score}")
    # ====================== Simulated Annealing ========================
    # for protein_sequence in protein_sequences:
    #     annealing = SimulatedAnnealing(protein_sequence)
    #     annealing.run(save_plot = True, save_data= True, repeats = 1)

    #     print(f"Best score for Simulated Annealing algorithm for protein {protein_sequence}: {annealing.best_score}")
    pass