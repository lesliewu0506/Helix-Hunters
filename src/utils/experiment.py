import src.visualisation as vis
from src.algorithms import Random, Greedy, HillClimber, SimulatedAnnealing
from src.utils import PROTEIN_SEQUENCES, ITERATIVE_ALGORITHM_FACTOR

algorithms = [("Random", Random), ("Greedy", Greedy), ("Hill Climber", HillClimber), ("Simulated Annealing", SimulatedAnnealing)]

def run(protein_sequence: str ="all", algorithm: str = "all", dimension: int = 3, show: bool = False, save: bool = True, repeats: int = 1, iterations: int = 10000) -> None:
    """
    Main function for running the experiment.
    It will run every algorithm on every protein sequence, save the data and print the best score.
    Or if a specific protein sequence is given, it will only run that protein sequence.
    Or if a specific algorithm is given, it will only run that algorithm.
    A combination of both is also possible.

    Parameters
    ----------
    algorithm: str, optional
        The algorithm to run. If no argument is given, all algorithm will be run. Default is `"all"`.

    dimension : int, optional
        The dimension in which the folding takes place (`2` or `3`). Default is `3`.

    show : bool, optional
            If `True` show the plot. Default is `False`.

    save : bool, optional
        If `True` save the plot and data. Default is `True`.

    repeats : int, optional
        The number of independent runs to perform. Default is `1`, due to run time.
        
    iterations : int, optional
        The number of iterations per run. Default is `10000`, due to run time.
    
    Notes
    -----
    `Hill Climber` and `Simulated Annealing` have default `iterations = 1000`.
    This is because their run time is extremely long for `10000` iterations.
    
    Raises
    ------
    ValueError
        - If an invalid algorithm name is given.
        - if an invalid protein sequence is given.
        - If an invalid dimension is given, must choose from `2` or `3`.
        - If an invalid repeats or iterations is given, must use >=1.

    Examples
    --------
    - Run all

    >>> run(dimension = 2, show = False, save = True, repeats = 3, iterations = 5000)

    - Run Greedy

    >>> run(algorithm = "Greedy", dimension = 3, repeats = 10, iterations = 1000)

    - Run "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"

    >>> run(protein_sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH", dimension = 3, show = True, save = False, repeats = 1, iterations = 10000)
    """
    valid_algorithms = ["Random", "Greedy", "Hill Climber", "Simulated Annealing", "all"]
    if algorithm not in valid_algorithms:
        raise ValueError(f"Invalid algorithm name. Choose from:\n{valid_algorithms}")
    
    valid_protein_sequences: list[str] = PROTEIN_SEQUENCES + ["all"]
    if protein_sequence not in valid_protein_sequences:
        raise ValueError(f"Invalid protein sequence. Choose from:\n{valid_protein_sequences}")

    if dimension not in [2, 3]:
        raise ValueError("Invalid dimension given. Choose from:\n[2, 3].")

    if repeats < 1 or iterations < 1:
        raise ValueError("Invalid repeats or iterations given. Use >= 1.")
    
    # Filter protein sequence
    selected_protein_sequence =  (PROTEIN_SEQUENCES if protein_sequence == "all" 
        else [seq for seq in PROTEIN_SEQUENCES if seq == protein_sequence])

    # Filter algorithm
    selected_algorithm = (algorithms if algorithm == "all" 
        else [alg for alg in algorithms if alg[0] == algorithm])
    
    for protein_sequence in selected_protein_sequence:
        for algorithm, cls in selected_algorithm:
            _run_algorithm(cls, algorithm, protein_sequence, dimension, show, save, repeats, iterations)

def _run_algorithm(
        cls: Random | Greedy | HillClimber | SimulatedAnnealing,
        algorithm: str,
        protein_sequence: str,
        dimension: int,
        show: bool,
        save: bool,
        repeats: int,
        iterations: int
        ) -> None:
    """Helper function for running an algorithm, saving the data and printing the best scores."""
    if cls in [SimulatedAnnealing, HillClimber]:
        iterations = iterations // ITERATIVE_ALGORITHM_FACTOR

    instance: Random | Greedy | HillClimber | SimulatedAnnealing = cls(protein_sequence, dimension)
    instance.run(show_plot = show, save_plot = save, save_data = save, repeats = repeats, iterations = iterations)
    print(f"Best score for {algorithm} algorithm for protein {protein_sequence}: {instance.best_score}")

def view(
        protein_sequence: str = "all",
        dimension: int = 3,
        show_plot: bool = True,
        save_plot: bool = True
        ) -> None:
    """
    Shows the boxplots for the different algorithms and saves the boxplots.
    Takes in protein sequence and dimension as arguments and plots the corresponding boxplot.
    If no arguments are given, it will plot all boxplots in 3D and show and save the plots.

    Parameters
    ----------
    protein_sequence : str, optional
        Protein sequence (for example `HHPHHHPH`). Default is `"all"`.

    dimension : int, optional
        The dimension in which the folding takes place (`2` or `3`). Default is `3`.

    show_plot : bool, optional
        If `True` show the plot. Default is `True`.

    save_plot : bool, optional
        If `True` save the plot. Default is `True`.

    Raises
    ------
    ValueError
        - If the protein sequence is invalid.
        - If an invalid dimension is given, must choose from `2` or `3`.

    Example
    -------
    >>> view(protein_sequence = "HHPHHHPH", dimension = 2, show_plot = True, save_plot = False)
    """
    # Check arguments given
    if protein_sequence != "all" and protein_sequence not in PROTEIN_SEQUENCES:
        raise ValueError(f"Invalid protein sequence. Choose from:\n{PROTEIN_SEQUENCES}")
    if dimension not in [2, 3]:
        raise ValueError("Invalid dimension given. Choose from:\n[2, 3].")

    if protein_sequence == "all":
        for sequence in PROTEIN_SEQUENCES:
            vis.boxplot(sequence, dimension, show_plot, save_plot)
    else:
        vis.boxplot(protein_sequence, dimension, show_plot, save_plot)