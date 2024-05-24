import numpy as np
import itertools

# Function to read coordinates from a file
def read_coordinates(filename):
    coordinates = []
    count = 0 
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.split())
            coordinates.append((x, y))
            count += 1
            if(count > 50):
                break
    return np.array(coordinates)

# Function to calculate the Euclidean distance matrix
def calculate_distance_matrix(coordinates):
    n = len(coordinates)
    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                dist_matrix[i][j] = np.linalg.norm(coordinates[i] - coordinates[j])
    return dist_matrix

# Branch and Bound class for TSP
class TSPSolver:
    def __init__(self, dist_matrix):
        self.n = len(dist_matrix)
        self.dist_matrix = dist_matrix
        self.final_res = float('inf')
        self.final_path = []

    # Function to copy temporary solution to final solution
    def copy_to_final(self, curr_path):
        self.final_path = curr_path[:]
        self.final_path.append(curr_path[0])

    # Function to find the minimum edge cost
    def first_min(self, i):
        min_val = float('inf')
        for k in range(self.n):
            if self.dist_matrix[i][k] < min_val and i != k:
                min_val = self.dist_matrix[i][k]
        return min_val

    # Function to find the second minimum edge cost
    def second_min(self, i):
        first, second = float('inf'), float('inf')
        for j in range(self.n):
            if i == j:
                continue
            if self.dist_matrix[i][j] <= first:
                second = first
                first = self.dist_matrix[i][j]
            elif self.dist_matrix[i][j] <= second and self.dist_matrix[i][j] != first:
                second = self.dist_matrix[i][j]
        return second

    # Function to solve TSP using Branch and Bound
    def tsp_branch_and_bound(self, curr_bound, curr_weight, level, curr_path, visited):
        if level == self.n:
            last_to_first = self.dist_matrix[curr_path[level - 1]][curr_path[0]]
            if last_to_first != 0:
                curr_res = curr_weight + last_to_first
                if curr_res < self.final_res:
                    self.copy_to_final(curr_path)
                    self.final_res = curr_res
            return

        for i in range(self.n):
            if (self.dist_matrix[curr_path[level-1]][i] != 0 and visited[i] == False):
                temp = curr_bound
                curr_weight += self.dist_matrix[curr_path[level - 1]][i]
                if level == 1:
                    curr_bound -= ((self.first_min(curr_path[level - 1]) + self.first_min(i)) / 2)
                else:
                    curr_bound -= ((self.second_min(curr_path[level - 1]) + self.first_min(i)) / 2)

                if curr_bound + curr_weight < self.final_res:
                    curr_path[level] = i
                    visited[i] = True
                    self.tsp_branch_and_bound(curr_bound, curr_weight, level + 1, curr_path, visited)

                curr_weight -= self.dist_matrix[curr_path[level - 1]][i]
                curr_bound = temp
                visited = [False] * self.n
                for j in range(level):
                    if curr_path[j] != -1:
                        visited[curr_path[j]] = True

    def solve(self):
        curr_bound = 0
        curr_path = [-1] * (self.n + 1)
        visited = [False] * self.n

        for i in range(self.n):
            curr_bound += (self.first_min(i) + self.second_min(i))

        curr_bound = np.ceil(curr_bound / 2)
        visited[0] = True
        curr_path[0] = 0

        self.tsp_branch_and_bound(curr_bound, 0, 1, curr_path, visited)

        return self.final_path, self.final_res

# Example usage
filename = 'points50_china.txt' 
coordinates = read_coordinates(filename)
dist_matrix = calculate_distance_matrix(coordinates)

solver = TSPSolver(dist_matrix)
path, cost = solver.solve()

print("Optimal Tour:", path)
print("Optimal Tour Cost:", cost)
