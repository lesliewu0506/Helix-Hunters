import click
import src.brute_force as brute

from src.utils import run, view, PROTEIN_SEQUENCES, ALGORITHMS

@click.group()
def cmd_group():
    """
    Main CLI for protein folding.

    \b
    1. run:
    Main function for running the experiment.
    It will run every algorithm on every protein sequence, save the data and print the best score.
    Or if a specific protein sequence is given, it will only run that protein sequence.
    Or if a specific algorithm is given, it will only run that algorithm.
    A combination of both is also possible.

    \b
    2. view:
    Main function for showing the boxplots for the different algorithms and saves the boxplots.
    Takes in protein sequence and dimension as arguments and plots the corresponding boxplot.
    If no arguments are given, it will plot all boxplots in 3D and show and save the plots.

    \b
    3. bruteforce:
    Main function for generating all protein foldings structures in 2D and evaluating all of them.
    Finds the protein structures with the lowest scores and saves them. This will only evaluate the
    first two protein sequences. 

    Usage:

    - Run all algorithms on all sequences in 3D:

        $ python main.py run
    
    - View results for all protein sequences in 3D:

        $ python main.py view

    - Brute force proteins sequences in 2D:

        $ python main.py bruteforce
    """
    pass

@cmd_group.command("run")
@click.option(
    "-p",
    "--protein",
    default = "all",
    show_default = True,
    type = str,
    help = f"Give the protein sequence to fold. Choose from:\n{PROTEIN_SEQUENCES}. If not given, all protein sequences will be ran."
)
@click.option(
    "-a",
    "--algorithm",
    default = "all",
    show_default = True,
    type = str,
    help = f"Give the name of the algorithm to use. Choose from:\n{ALGORITHMS}. If not given, all algorithms will be ran."
)
@click.option(
    "--graph",
    is_flag = True,
    show_default = True,
    help = "Shows the plots."
)
@click.option(
    "--save",
    is_flag = True,
    show_default = True,
    help = "Saves the plots and data generated."
)
@click.option(
    "-d",
    "--dimension",
    default = 3,
    show_default = True,
    type = click.IntRange(2, 3),
    help = "Dimension in which the folding takes place (either 2 or 3)."
)
@click.option(
    "-r",
    "--repeats",
    default = 1,
    show_default = True,
    type = click.IntRange(1, 20),
    help = "Number of repeats (must be > 0) per run."
)
@click.option(
    "-i",
    "--iterations",
    default = 10000,
    show_default = True,
    type = int,
    help = "Number of iterations per repeat (must be > 0)."
)
def run_experiment(
    protein: str,
    algorithm: str,
    dimension: int,
    graph: bool,
    save: bool,
    repeats: int,
    iterations: int) -> None:
    """
    Main function for running the experiment.

    \b
    It will run every algorithm on every protein sequence, save the data and print the best score.
    Or if a specific protein sequence is given, it will only run that protein sequence.
    Or if a specific algorithm is given, it will only run that algorithm.
    A combination of both is also possible.
    
    Notes

    Hill Climber and Simulated Annealing have default iterations = 1000.
    This is because their run time is extremely long for 10000 iterations.
    
    Examples:

    1. Run all algorithms on all protein sequences in 3D:

        $ python main.py run

    2. Run Hill Climber on a specific protein sequence with 5 repeats in 2D, showing the plot and 1000 iterations:

    \b
        $ python main.py run -p HHPHHHPHPHHHPH -a "Hill Climber" -d 2 --graph -r 5 -i 1000
    """
    run(
        protein_sequence = protein,
        algorithm = algorithm,
        dimension = dimension,
        show = graph,
        save = save,
        repeats = repeats,
        iterations = iterations)
@cmd_group.command("view")
@click.option(
    "-p",
    "--protein",
    default = "all",
    show_default = True,
    type = str,
    help = f"Give the protein sequence to fold. Choose from:    {PROTEIN_SEQUENCES}. If not given, all protein sequences will be ran."
)
@click.option(
    "--graph",
    is_flag = True,
    show_default = True,
    help = "Shows the plots."
)
@click.option(
    "--save",
    is_flag = True,
    show_default = True,
    help = "Saves the plots and data generated."
)
@click.option(
    "-d",
    "--dimension",
    default = 3,
    show_default = True,
    type = click.IntRange(2, 3),
    help = "Dimension in which the folding takes place (either 2 or 3)."
)
def view_experiment(
    protein: str,
    dimension: int,
    graph: bool,
    save = bool) -> None:
    """
    Main function for viewing the results of the experiment.

    \b
    Shows the boxplots for the different algorithms and saves the boxplots.
    Takes in protein sequence and dimension as arguments and plots the corresponding boxplot.
    If no arguments are given, it will plot all boxplots in 3D and show and save the plots.
    
    Examples:

    1. Show all boxplots for all protein sequences:

        $ python main.py view

    2. Show boxplot for specific protein sequence in 2D, show and save it:

    \b
        $ python main.py view -p HPHPPHHPHPPHPHHPPHPH -d 2 --graph --save
    
    """
    view(
        protein_sequence = protein,
        dimension = dimension,
        show_plot = graph,
        save_plot = save)
    
@cmd_group.command("bruteforce")
def brute_force():
    """
    Main function for brute forcing all possible folding structures.
    
    It generates all possible folding combinations into a CSV file.
    It then reads all the structures from a CSV file.
    Using the multiprocessing module, it evaluates every structure with helper functions.
    Plots and saves the best structure and prints the rating of the best structure.
    Finally saves the directions into a CSV file.
    """
    for protein_sequence in PROTEIN_SEQUENCES[:2]:
        brute.generate_all_foldings(protein_sequence)
        brute.brute_force(protein_sequence, save = True)