import random as rd
import numpy as np
from typing import List, Callable, Tuple

class Genetic_Algorithm():
    """The Genetic random class generates a sequence for the folding direction."""

    def __init__(self, protein_sequence: str, population_size: int, mutation_rate: float) -> None:
        self.protein_sequence = protein_sequence
        self.population_size = population_size
        self.mutation_rate = mutation_rate

    def run(self, show_plot: bool = False, save_plot: bool = False, save_data: bool = False, repeats: int = 1, iterations: int = 10000) -> None:
        """The Genetic algorithm generates sequences for a protein and calculates the scores
        self.run_algorithm(
            algorithm = "Genetic",
            show_plot = show_plot,
            save_plot = save_plot,
            saved_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._greedy_iterated)
        
    def initialized_population(self):
        population = list()
        target_length = len(target)

        for i in range(population_size):
            temperature = list()
            for j in range(target_length):
                temperature.append(random.choice(Genes))
            population.append(temperature)

        return population
    
    def generate_crossover(selected_chromo, chromo_length, population):
        offspring_cross = []
        for i in range(int(population_size)):
            parent_1= random.choice(selected_chromo)
            parent_2 = random.choice(population[:int(population_size*50)])

            p_1 = parent_1[0]
            p_2 = parent_2[0]

            cross_point = random.randit(1, chromo_length-1)
            child = p_1[:crossover_point] + p_2[crossover_point:]
            offspring_cross.extend([child])
        return offspring_cross
    
    def apply_mutate(offspring, mutation_rate):
        mutated_offspring = []



"""

def create_initial_population(size, lower_bound, upper_bound):


initialized_population()
evaluate_fit()
select_parents()
generate_crossover()
apply_mutation()
check_convergence()