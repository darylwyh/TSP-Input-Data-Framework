from problem_tsp import TSP
import numpy as np
import random
import matplotlib.pyplot as plt
from tsp_baseline import run_insertion

def load_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.split())
            points.append((x, y))
    return np.array(points)

def main():
    # Load points from the file
    filename = 'attention-learn-to-route\data\TSPDataset\points20_second.txt'
    points = load_points(filename)

    # Run Random Insertion algorithm
    cost, tour = run_insertion(points, method='random')

    # Print the results
    print(f"Total cost: {cost:.2f}")
    print(f"Tour: {tour}")

if __name__ == "__main__":
    main()