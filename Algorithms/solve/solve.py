import numpy as np
from tsp_baseline import run_insertion, solve_lkh_log, load_dataset, save_dataset, get_lkh_executable

def load_points(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.split())
            points.append((x, y))
    return np.array(points)

def main():
    # n=10,20,30,40,50 
    # Load points from the file
    filename = 'points20_second.txt'
    points = load_points(filename)

    # Run Random Insertion algorithm
    cost, tour = run_insertion(points, method='random')
    # Run 
    points = 10
    # Print the results
    print(f"N is: {points}")
    print(f"Total cost: {cost:.2f}")
    print(f"Tour: {tour}")

if __name__ == "__main__":
    main()
