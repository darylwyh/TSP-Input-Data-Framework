import numpy as np
from geopy.distance import geodesic
 
from sklearn.manifold import MDS


# Define the file path
file_path = 'city_coordinates.txt'
# Read city coordinates from the file
city_coordinates = {}
with open(file_path, 'r') as file:
    for line in file:
        parts = line.strip().split(',')
        city_name = parts[0]
        latitude = float(parts[1])
        longitude = float(parts[2])
        city_coordinates[city_name] = (latitude, longitude)

# Function to parse the line from the weights.intra file
def parse_line(line):
    parts = line.split()
    city1 = parts[0].split(',')[0]
    city2 = parts[1].split(',')[0]
    weight = float(parts[2])  # Assuming weight is a placeholder for distance here
    return city1, city2, weight

# Read city pairs from file and create dictionary
city_pairs = {}
with open('weights.intra', 'r') as file:
    input_data = file.read().splitlines()

for line in input_data:
    city1, city2, weight = parse_line(line)
    city_pairs[(city1, city2)] = weight

# Create a list of unique cities
cities = list(city_coordinates.keys())

# Create a distance matrix initialized with zeros
num_cities = len(cities)
distance_matrix = np.zeros((num_cities, num_cities))

# Compute geodesic distances and fill the distance matrix
for i, city1 in enumerate(cities):
    for j, city2 in enumerate(cities):
        if i != j:
            # Calculate real-world geodesic distance
            distance_matrix[i][j] = geodesic(
                city_coordinates[city1],
                city_coordinates[city2]
            ).kilometers

# Normalize the distance matrix
max_distance = np.max(distance_matrix)
distance_matrix = distance_matrix / max_distance

print("Normalized Distance Matrix:")
print(distance_matrix)

# where 'distance_matrix' is your (N x N) distance matrix
embedding = MDS(n_components=2, dissimilarity='precomputed')
points = embedding.fit_transform(distance_matrix)
points = points[:20]  # Select only 100 points if there are more
# assert points.shape == (100, 2), "The shape of the points array must be (100, 2)"

print(points)

# Save to file
np.savetxt('points20.txt', points, fmt='%f')
'''
# Simulated city coordinates mapping; replace with your actual data
city_coordinates = {
    'Darwin': (-12.4634, 130.8456),
    'Perth': (-31.9505, 115.8605),
    'Sydney': (-33.8688, 151.2093),
    'Brisbane': (-27.4698, 153.0251),
    'Canberra': (-35.2809, 149.1300)
}

'''