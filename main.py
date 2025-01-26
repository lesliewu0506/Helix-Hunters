# Import Brute Force Scripts
from src.brute_force import brute_force, generate_all_foldings

# Import Helpers 
from src.utils import run, view, PROTEIN_SEQUENCES

if __name__ == "__main__":
    # =========================== Experiment ===========================
    # This will collect data for different algorithms
    run(protein_sequence = "all", algorithm = "all", show = False, save = True, dimension = 3, repeats = 20, iterations = 10000)

    # =========================== Visualisation ========================
    # View the boxplots for the different distributions
    view(protein_sequence = "all", dimension = 3, show_plot = True, save_plot = True)

    # =========================== Brute Force ==========================
    # NOTE: This function is good only for the first two protein sequences in 2D.
    # Running it for the others results in huge data files (>>5GB). 
    # This is only to showcase how hard it is to try all possible states.

    # Generate all possible folding structures
    # for protein_sequence in PROTEIN_SEQUENCES[:2]:
    #     generate_all_foldings(protein_sequence)
    #     # Try all possible combinations
    #     brute_force(protein_sequence, save = True)