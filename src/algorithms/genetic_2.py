import random as rd
import matplotlib as plt
import numpy as np
from prettytable import PrettyTable

from src.utils.helpers import random_fold
from src.classes.protein import Protein
from typing import Callable, Optional
from .general import General

class Genetic_Algorithm():
    """The Genetic random class generates a sequence for the folding direction."""

    def __init__(self, protein_sequence: str) -> None:
        super().__init__(protein_sequence)

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
        
    def initialize_population(target):
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

def evaluated_fit(parameters):
    a, b, c = parameters
    if a <= 0:
        return = -float('inf')
    vertex_x = -b / (2*a) 
    vertex_y = a * (vertex_x ** 2) + b * vertex_x + c
    y_left = a * (-1) ** 2 + b * (-1) + c
    y_right = a * (-1) ** 2 + b * (1) + c
    curviness = abs(y_left - vertex_y) + abs(y_right - vertex_y)
    return -curviness

def create_initial_population(size, lower_bound, upper_bound):


initialized_population()
evaluate_fit()
select_parents()
generate_crossover()
apply_mutation()
check_convergence()