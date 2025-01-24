# Import Brute Forcing Scripts
from src.brute_force import brute_force, generate_all_foldings

# Import Algorithms
from src.algorithms import Random, Greedy, HillClimber, SimulatedAnnealing

# Import Helpers 
from src.utils import run, view, protein_sequences

if __name__ == "__main__":
    # =========================== Brute Force ==========================
    # NOTE: This function is good only for the first three protein sequences.
    # Running it for the others results in huge data files (>5GB). 
    # This is only for demonstration purposes.

    # Generate all possible folding structures
    # for protein_sequence in protein_sequences[:3]:
    #     generate_all_foldings(protein_sequence)
    #     # Try all possible combinations
    #     brute_force(protein_sequence, save = True)

    # ==================================================================
    # Experiment
    # ==================================================================
    # This will collect data for different algorithms.
    # run(dimension = 3, repeats = 1, iterations = 10000)
    # View the boxplots for the different distributions
    # view(dimension = 2, protein_sequence = "all")

    # ==================================================================
    # The following sections generate data for each algorithm seperately

    # # =========================== Random ===============================
    # for protein_sequence in protein_sequences:
    #     random = Random(protein_sequence, dimension = 3)
    #     random.run(show_plot = False, save_plot = True, save_data = True, repeats = 1, iterations = 10000)

    #     print(f"Best score for Random algorithm for protein {protein_sequence}: {random.best_score}")

    # # ====================== Random Greedy ==============================
    # for protein_sequence in protein_sequences:
    #     greedy = Greedy(protein_sequence, dimension = 3)
    #     greedy.run(show_plot = False, save_plot = True, save_data = True, repeats = 1, iterations = 10000)

    #     print(f"Best score for Greedy algorithm for protein {protein_sequence}: {greedy.best_score}")

    # # ====================== Hill Climber ===============================
    # for protein_sequence in protein_sequences:
    #     hillclimber = HillClimber(protein_sequence, dimension = 3)
    #     hillclimber.run(show_plot = False, save_plot = True, save_data= True, repeats = 1, iterations = 1000)

    #     print(f"Best score for Hill Climber algorithm for protein {protein_sequence}: {hillclimber.best_score}")

    # # ====================== Simulated Annealing ========================
    # for protein_sequence in protein_sequences:
    #     annealing = SimulatedAnnealing(protein_sequence, dimension = 3)
    #     annealing.run(show_plot = False, save_plot = True, save_data= True, repeats = 1, iterations = 1000)

    #     print(f"Best score for Simulated Annealing algorithm for protein {protein_sequence}: {annealing.best_score}")
    random = Random(protein_sequences[-1], dimension = 3)
    random.run(show_plot = True, iterations = 1000)
    pass