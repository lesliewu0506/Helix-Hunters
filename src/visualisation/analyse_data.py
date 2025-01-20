import csv
import matplotlib.pyplot as plt

def import_data(protein_sequence: str, algorithm: str, folder: str):
    histogram_data: list[list[int]] = []

    with open(f"data/histogram_data/{folder}/{algorithm}_{protein_sequence}.csv", "r") as csvfile:

        reader = csv.reader(csvfile)
        
        for row in reader:
            histogram: list[int] = [int(value) for value in row]
            histogram_data.append(histogram)

    csvfile.close()

    return histogram_data
