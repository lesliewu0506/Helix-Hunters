from src.algorithms.randomise import Random 
from src.algorithms.greedy import Greedy
from src.algorithms.hill_climber import HillClimber
from src.algorithms.simulated_annealing import SimulatedAnnealing
from src.visualisation.analyse_data import boxplot
from src.utils.helpers import protein_sequence_map

protein_sequences = ["HHPHHHPHPHHHPH",
                     "HPHPPHHPHPPHPHHPPHPH",
                     "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
                     "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
                     "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
                     "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
                     "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
                     "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]

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

    
def view():
    for protein_sequence in protein_sequences:
        boxplot(protein_sequence, protein_sequence_map[protein_sequence])