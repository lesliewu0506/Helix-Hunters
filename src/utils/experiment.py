from src.algorithms import Random, Greedy, HillClimber, SimulatedAnnealing
from src.visualisation import boxplot
from src.utils.constants import protein_sequence_map, protein_sequences

def run(repeats: int = 1, iterations: int = 10000) -> None:
    """Helper function for collecting data for all algorithms."""
    for protein_sequence in protein_sequences:

        random = Random(protein_sequence)
        random.run(save_plot = True, save_data = True, repeats = repeats, iterations = iterations)
        print(f"Best score for Random algorithm for protein {protein_sequence}: {random.best_score}")

        greedy = Greedy(protein_sequence)
        greedy.run(save_plot = True, save_data = True, repeats = repeats, iterations = iterations)
        print(f"Best score for Greedy algorithm for protein {protein_sequence}: {greedy.best_score}")

        hillclimber = HillClimber(protein_sequence)
        hillclimber.run(save_plot = True, save_data= True, repeats = repeats, iterations = iterations / 10)
        print(f"Best score for Hill Climber algorithm for protein {protein_sequence}: {hillclimber.best_score}")

        annealing = SimulatedAnnealing(protein_sequence)
        annealing.run(save_plot = True, save_data= True, repeats = repeats, iterations = iterations / 10)
        print(f"Best score for Simulated Annealing algorithm for protein {protein_sequence}: {annealing.best_score}")

    
def view(protein_sequence: str = "all") -> None:
    """
    Shows the boxplots for the different algorithms and saves the boxplots.
    Takes in a protein sequence as argument and plots the corresponding boxplot.
    If no argument is given, plots all boxplots.
    """
    if protein_sequence == "all":
        for sequence in protein_sequences:
            boxplot(sequence, protein_sequence_map[sequence])
    else:
        boxplot(protein_sequence, protein_sequence_map[protein_sequence])