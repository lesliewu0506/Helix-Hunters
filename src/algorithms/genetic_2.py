import os
import random as rd
import csv
from typing import List, Optional

class Protein:
    """This class generates a sequence for the folding direction."""
    
    def __init__(self, protein_sequence: str, amino_directions: Optional[List[int]] = None) -> None:
        self.protein_sequence: str = protein_sequence
        self.amino_directions: Optional[List[int]] = amino_directions
        self.structure: Optional[dict] = None
        self.protein_rating: int = 1

    def build_no_function(self) -> None:
        """Builds structure"""
        # Placeholder implementation for building a structure
        if self.amino_directions:
            self.structure = {i: (amino, direction) for i, (amino, direction) in enumerate(zip(self.protein_sequence, self.amino_directions))}

    def get_rating(self) -> int:
        """Return the rating of the protein fold."""
        return self.protein_rating
    
class Genetic_Algorithm:

    def __init__(self, protein_sequence: str):
        self.protein_sequence: str = protein_sequence
        self.histogram_data: List[List[int]] = []
        self.best_score: int = 0
        self.best_protein: Optional[Protein] = None

        self.protein_sequence_map = {
            "HHPHHHPHPHHHPH": "1",
            "HPHPPHHPHPPHPHHPPHPH": "2",
            "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP": "3",
            "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH": "4",
            "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP": "5",
            "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC": "6",
            "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH": "7",
            "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH": "8"
        }
        self.folder = self.protein_sequence_map.get(protein_sequence, "default")
        self.output_dir = f"c:/mnt/Users/trins/Downloads/{self.folder}/genetic_{self.protein_sequence}.csv"
        os.makedirs(self.output_dir, exist_ok=True)

    def run(self, iterations: int = 10000, population_size: int = 100, mutation_rate: float = 0.01) -> None:
        """Run the Genetic Algorithm for 10000 iterations."""
        population = self.Initialized_population(population_size)

        for iteration in range(iterations):
            fitness_scores = self.evaluated_population(population)
            best_index = fitness_scores.index(max(fitness_scores))
            best_individual = population[best_index]

            if self.best_protein is None or self.evaluated_fitness(best_individual) > self.best_score:
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

    def Initialized_population(self, population_size: int) -> List[List[int]]:
        """Initialize the population with random sequences."""
        population = []
        target_length = len(self.protein_sequence)
        genes = [0, 1, 2, 3]

        for _ in range(population_size):
            individual = [rd.choice(genes) for _ in range(target_length)]
            population.append(individual)

        return population

    def evaluated_population(self, population: List[List[int]]) -> List[float]:
        """Evaluating the fitness of each individual in the population."""
        return [self.evaluated_fitness(individual) for individual in population]

    def evaluated_fitness(self, individual: List[int]) -> float:
        """Evaluating the fitness of an individual sequence."""
        protein = Protein(self.protein_sequence, individual)
        protein.build_no_function()
        return protein.get_rating()

    def select_parents(self, population: List[List[int]], fitness_scores: List[float]) -> List[List[int]]:
        """Select parents based on fitness scores using the selection method."""
        parents = []
        for _ in range(len(population)):
            candidates = rd.sample(list(zip(population, fitness_scores)), 3)
            best_candidate = max(candidates, key=lambda x: x[1])
            parents.append(best_candidate[0])
        return parents
    
    def generate_crossover(self, parents: List[List[int]]) -> List[List[int]]:
        """Generate offspring using crossover."""
        offspring = []
        for _ in range(len(parents) // 2):
            parent1, parent2 = rd.sample(parents, 2)
            crossover_point = rd.randint(1, len(parent1) - 1)
            child1 = parent1[:crossover_point] + parent2[crossover_point:]
            child2 = parent2[:crossover_point] + parent1[crossover_point:]
            offspring.extend([child1, child2])
        return offspring

    def apply_mutation(self, offspring: List[List[int]], mutation_rate: float) -> List[List[int]]:
        """Apply mutation to offspring."""
        mutated_offspring = []
        for individual in offspring:
            mutated_individual = [
                gene if rd.random() > mutation_rate else rd.choice([0, 1, 2, 3])
                for gene in individual
            ]
            mutated_offspring.append(mutated_individual)
        return mutated_offspring

    def check_convergence(self, fitness_scores: List[float]) -> bool:
        """Check if the population convergence has been acquired."""
        return len(set(fitness_scores)) == 1

    def output_csv(self) -> None:
        """Saves histogram data into a csv file."""
        csv_path = os.path.join(self.output_dir, f"genetic_{self.protein_sequence}.csv")
        os.makedirs(self.output_dir, exist_ok=True) 

        with open(csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            for histogram in self.histogram_data:
                writer.writerow(histogram)

        print(f"Histogram data saved to: {csv_path}")

# Output
if __name__ == "__main__":
    ga = Genetic_Algorithm(protein_sequence="HHPHHHPHPHHHPH")
    ga.run(iterations=1000, population_size=200, mutation_rate=0.05)
    ga.output_csv()