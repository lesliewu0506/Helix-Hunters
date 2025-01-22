# Import Brute Forcing Scripts
from src.brute_force import brute_force, generate_all_foldings

# Import Algorithms
from src.algorithms import Random, Greedy, HillClimber, SimulatedAnnealing

# Import Helpers 
from src.utils import run, view
from src.utils.constants import protein_sequences

if __name__ == "__main__":
    # =========================== Brute Force ==========================
    # NOTE: This function is good only for the first protein sequence.
    # Running it for the others results in huge data files (>5GB). 
    # This is only for demonstration purposes.

    # # Generate all possible folding structures
    # generate_all_foldings(protein_sequences[0])
    # # Try all possible combinations
    # brute_force(protein_sequences[0], save = True)

    # ==================================================================
    # Experiment
    # ==================================================================
    # This will collect data for different algorithms.
    # # Repeats is how many times one algorithm should run.
    # # Per run there are then iterations amount of iterations.
    # run(repeats = 10, iterations = 10000)
    # # View the boxplots for the different distributions
    # view(protein_sequences[-1])

    # ==================================================================
    # The following functions generate data for each algorithm seperately

    # =========================== Random ===============================
    # for protein_sequence in protein_sequences:
    #     random = Random(protein_sequence)
    #     random.run(save_plot = True, save_data = True, repeats = 10)

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
        # annealing = SimulatedAnnealing(protein_sequence)
        # annealing.run(save_plot = True, save_data= True, repeats = 1)

    #     print(f"Best score for Simulated Annealing algorithm for protein {protein_sequence}: {annealing.best_score}")
    random = Random(protein_sequences[0], 2)
    random.run(show_plot = True)
    # annealing = SimulatedAnnealing(protein_sequences[0], 2)
    # annealing.run(show_plot=True, repeats = 1, iterations = 1)
    # hill = HillClimber(protein_sequences[0], 3)
    # hill.run(show_plot=True)
    # print(hill.best_protein.amino_directions)
    pass