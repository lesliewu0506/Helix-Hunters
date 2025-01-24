# Import Brute Force Scripts
from src.brute_force import brute_force, generate_all_foldings

# Import Helpers 
from src.utils import run, view, protein_sequences

if __name__ == "__main__":
    # =========================== Experiment ===========================
    # This will collect data for different algorithms
    run(protein_sequence = "all", algorithm = "all", show = False, save = False, dimension = 3, repeats = 1, iterations = 1000)
    # View the boxplots for the different distributions
    view(protein_sequence = "all", dimension = 3, show_plot = True, save_plot = False)

    # =========================== Brute Force ==========================
    # NOTE: This function is good only for the first three protein sequences.
    # Running it for the others results in huge data files (>>5GB). 
    # This is only for demonstration purposes.

    # Generate all possible folding structures
    # for protein_sequence in protein_sequences[:2]:
    #     generate_all_foldings(protein_sequence)
    #     # Try all possible combinations
    #     brute_force(protein_sequence, save = True)