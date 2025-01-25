import random as rd

from . import General
from src.classes import Protein
from src.utils import random_fold, POPULATION_SIZE, MUTATION_RATE, DIRECTION_CHOICES_3D

class Genetic(General):
    """
    The Genetic Algorithm class implements the optimization technique inspired by natural selection.
    It generates solutions by combining and mutating parent solutions.
    The fittest solutions are selected for the next generation.
    """

    def __init__(
            self,
            protein_sequence: str,
            dimension: int,
            population_size: int = POPULATION_SIZE,
            mutation_rate: float = MUTATION_RATE
            ) -> None:
        if dimension not in [2, 3]:
            raise ValueError("Invalid dimension given. Choose from:\n[2, 3].")

        super().__init__(protein_sequence, dimension)
        self.population_size: int = population_size
        self.mutation_rate: float = mutation_rate

    def run(
        self,
        show_plot: bool = False,
        save_plot: bool = False,
        save_data: bool = False,
        repeats: int = 1,
        iterations: int = 10000
        ) -> None:
        """
        Runs the Genetic algorithm to optimize the protein structure.

        Parameters
        ----------
        show_plot : bool, optional
            If `True` show the plot. Default is `False`.

        save_plot : bool, optional
            If `True` save the plot. Default is `False`.

        save_data : bool, optional
            If `True`, saves the optimization results to a file. Default is `False`.

        repeats : int, optional
            The number of independent runs to perform. Default is `1`.

        iterations : int, optional
            The number of iterations per run. Default is `10000`.
        
        Raises
        ------
        ValueError
            If repeats or iterations has an invalid value (<1).
        """
        if repeats < 1 or iterations < 1:
            raise ValueError("Both repeats and iterations must be at least 1.")
        
        self.run_algorithm(
            algorithm = "Genetic",
            show_plot = show_plot,
            save_plot = save_plot,
            save_data = save_data,
            repeats = repeats,
            iterations = iterations,
            algorithm_function = self._genetic,
            population_size = self.population_size,
            mutation_rate = self.mutation_rate)

    def _genetic(self, iterations: int, population_size: int, mutation_rate: float, accept_function: None, temperature: float) -> None:
        """Run the Genetic Algorithm."""
        population: list[list[int]] = self.Initialized_population(population_size)

        for iteration in range(iterations):
            fitness_scores: list[int] = self.evaluated_population(population)
            best_index: int = fitness_scores.index(min(fitness_scores))
            best_individual: list[int] = population[best_index]

            if self.best_protein is None or self.evaluated_fitness(best_individual) < self.best_score:
                self.best_protein = Protein(self.protein_sequence, best_individual)
                self.best_protein.build_no_function()
                self.best_score = self.best_protein.get_rating()

            parents = self.select_parents(population, fitness_scores)
            offspring = self.generate_crossover(parents)
            population = self.apply_mutation(offspring, mutation_rate)

            self.histogram_data.append(fitness_scores)
            print(f"Iteration {iteration}: Best Score = {self.best_score}")

            if self.check_convergence(fitness_scores):
                print(f"Converged at iteration {iteration}")
                break

        print("Best protein fold:")
        print(self.best_protein.amino_directions)
        print(f"Best rating: {self.best_score}")

    def Initialized_population(self, population_size: int) -> list[list[int]]:
        """Initialize the population with random sequences."""
        population: list[list[int]] = []
        for _ in range(population_size):
            directions: list[int] = random_fold(self.protein_sequence, self.dimension)
            individual: Protein = Protein(self.protein_sequence, directions)
            individual.build_no_function()
            while individual.get_rating() == 1:
                directions = random_fold(self.protein_sequence, self.dimension)
                individual = Protein(self.protein_sequence, directions)
                individual.build_no_function()
            population.append(directions)

        return population

    def evaluated_population(self, population: list[list[int]]) -> list[float]:
        """Evaluating the fitness of each individual in the population."""
        return [self.evaluated_fitness(individual) for individual in population]

    def evaluated_fitness(self, individual: list[int]) -> float:
        """Evaluating the fitness of an individual sequence."""
        protein = Protein(self.protein_sequence, individual)
        protein.build_no_function()
        return protein.get_rating()

    def select_parents(self, population: list[list[int]], fitness_scores: list[float]) -> list[list[int]]:
        """Select parents based on fitness scores using the selection method."""
        parents = []
        for _ in range(len(population)):
            candidates = rd.sample(list(zip(population, fitness_scores)), 3)
            best_candidate = max(candidates, key=lambda x: x[1])
            parents.append(best_candidate[0])

        return parents
    
    def generate_crossover(self, parents: list[list[int]]) -> list[list[int]]:
        """Generate offspring using crossover."""
        offspring = []
        for _ in range(len(parents) // 2):
            parent1, parent2 = rd.sample(parents, 2)
            crossover_point = rd.randint(1, len(parent1) - 1)
            child1 = parent1[:crossover_point] + parent2[crossover_point:]
            child2 = parent2[:crossover_point] + parent1[crossover_point:]
            offspring.extend([child1, child2])
        return offspring

    def apply_mutation(self, offspring: list[list[int]], mutation_rate: float) -> list[list[int]]:
        """Apply mutation to offspring."""
        mutated_offspring = []
        for individual in offspring:
            mutated_individual = [
                gene if rd.random() > mutation_rate else rd.choice(DIRECTION_CHOICES_3D)
                for gene in individual
            ]
            mutated_offspring.append(mutated_individual)
        return mutated_offspring

    def check_convergence(self, fitness_scores: list[float]) -> bool:
        """Check if the population convergence has been acquired."""
        return len(set(fitness_scores)) == 1