import csv
import random
import numpy as np
from src.classes.grid import Grid
from src.classes.rating import Rating
from typing import List, Callable, Optional

class Genetic_Algorithm():
    """The Genetic random class generates a sequence for the folding direction."""

    def __init__(self, protein_sequence: str, population_size: int, mutation_rate: float) -> None:
        self.protein_sequence: str = protein_sequence
        self.population_size: Optional[List[int]] = amino_directions
        self.structure: Optional[Grid] = None
        self.protein_rating: int = 1

    def run(self, interations: int = 10000) -> None:
        """Run the Algorithm for iterations"""
        population = self.initialized_population()
        for iteration in range(interations):
            fitness_scores = self.evaluate_population(population)
            parents = self.select_parents(population, fitness_scores)
            offspring = self.generate_crossover(parents)
            population = self.apply_mutation(offspring)
            if self.check_convergence(fitness_scores):
                print(f"Convergence acquired at iteration {iteration}")
                break 
    
    def structure(self, function: Callable[[str], List[int]]) -> None:
        self.amino_directions = function(self.protein_sequence)
        structure = Grid(self.protein_sequence, self.amino_directions)

        # If structure not valid rating 1 
        if not structure.create_structure():
            self.protein_rating = 1
        else:
            self.structure = structure
            self.protein_rating = Rating(self.structure.get_structure()).get_rating()

    def building_structure(self) -> None:
        if self.amino_directions is not None:
            structure = Grid(self.protein_sequence, self.amino_directions)

            # Check if structure is valid else give rating 1
            if not structure.create_structure():
                self.protein_rating = 1
            else:
                self.structure = structure
                self.protein_rating = Rating(self.structure.get_structure()).get_rating()

       def output_csv(self, file_path: str = "output") -> None:
        """
        Creates a CSV file withamino acids and their fold.
        """
        if self.amino_directions is not None:
            with open(f'{file_path}.csv', 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)

                writer.writerow(['amino', 'fold'])

                for amino, direction in zip(self.protein_sequence, self.amino_directions):
                    writer.writerow([amino, direction])

                writer.writerow(['score', self.protein_rating])

    def initialized_population(self) -> List[List[int]]:
        """Initalize the population with sequences"""
        population = []
        target_length = len(self.population_size):
        genes = [0,1,2,3]

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



def create_initial_population(size, lower_bound, upper_bound):


initialized_population()
evaluate_fit()
select_parents()
generate_crossover()
apply_mutation()
check_convergence()