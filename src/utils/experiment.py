from src.algorithms import Random, Greedy, HillClimber, SimulatedAnnealing
from src.visualisation.boxplot import boxplot
from src.utils.constants import protein_sequences

def run(dimension: int = 2, repeats: int = 1, iterations: int = 10000) -> None:
    """
    Main function for running the experiment.
    It will run every algorithm on every protein sequence, save the data and print the best score.
    
    Parameters
    ----------
    dimension : int, optional
        The dimension in which the folding takes place (`2` or `3`). Default is `3`.

    repeats : int, optional
        The number of independent runs to perform. Default is `1`, due to run time.
        
    iterations : int, optional
        The number of iterations per run. Default is `1000`, due to run time.
    
    Notes
    -----
    `Hill Climber` and `Simulated Annealing` have default `iterations = 1000`.
    This is because their run time is extremely long for `10000` iterations.
    """
    for protein_sequence in protein_sequences:

        random = Random(protein_sequence, dimension)
        random.run(save_plot = True, save_data = True, repeats = repeats, iterations = iterations)
        print(f"Best score for Random algorithm for protein {protein_sequence}: {random.best_score}")

        greedy = Greedy(protein_sequence, dimension)
        greedy.run(save_plot = True, save_data = True, repeats = repeats, iterations = iterations)
        print(f"Best score for Greedy algorithm for protein {protein_sequence}: {greedy.best_score}")

        hillclimber = HillClimber(protein_sequence, dimension)
        hillclimber.run(save_plot = True, save_data= True, repeats = repeats, iterations = iterations // 10)
        print(f"Best score for Hill Climber algorithm for protein {protein_sequence}: {hillclimber.best_score}")

        annealing = SimulatedAnnealing(protein_sequence, dimension)
        annealing.run(save_plot = True, save_data= True, repeats = repeats, iterations = iterations // 10)
        print(f"Best score for Simulated Annealing algorithm for protein {protein_sequence}: {annealing.best_score}")

    
def view(protein_sequence: str = "all") -> None:
    """
    Shows the boxplots for the different algorithms and saves the boxplots.
    Takes in a protein sequence as argument and plots the corresponding boxplot.
    If no argument is given, it will plot all boxplots.

    Parameter
    ---------
    protein_sequence : str, optional
        Protein sequence (for example `HHHPPPHPCCP`). Default is `"all"`
    """
    if protein_sequence == "all":
        for sequence in protein_sequences:
            boxplot(sequence)
    else:
        boxplot(protein_sequence)