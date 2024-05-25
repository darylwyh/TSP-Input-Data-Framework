import pickle
import numpy as np

# Provided data
data = """
0.835065 0.490891
0.999999 0.273785
0.409980 0.670495
0.616275 0.076578
0.011318 0.155301
0.000001 0.137583
0.065130 0.074024
0.785564 0.000279
0.616868 0.060052
0.756233 0.006938
0.626757 0.088771
0.933133 0.381941
0.908322 0.116459
0.951923 0.104956
0.938414 0.107515
0.767353 0.028207
0.674327 0.000001
0.966952 0.129477
0.955299 0.116351
0.895653 0.066812
"""

# Convert the data to a list of tuples
points = [tuple(map(float, line.split())) for line in data.strip().split('\n')]

# Convert to a NumPy array
points_array = np.array(points)

# Save the data to a .pkl file
with open('points20_second.pkl', 'wb') as f:
    pickle.dump(points_array, f, pickle.HIGHEST_PROTOCOL)

print("Data has been saved to points20.pkl")
