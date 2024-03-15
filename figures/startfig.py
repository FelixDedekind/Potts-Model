import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Number of grid lines
num_lines = 4

# Create grid coordinates
x = np.linspace(0, 1, num_lines)
y = np.linspace(0, 1, num_lines)
z = np.zeros(num_lines)

# Generate random directions
directions = np.random.randint(0, 4, size=(num_lines, num_lines))

# Coordinates of vertices relative to the center of mass
vertices_relative = np.array([
    [1, 0, -1/np.sqrt(2)],
    [-1, 0, -1/np.sqrt(2)],
    [0, 1, 1/np.sqrt(2)],
    [0, -1, 1/np.sqrt(2)]
])

# Create figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot grid lines
for i in range(num_lines):
    ax.plot(x, np.ones(num_lines) * y[i], zs=z[i], color='black', alpha=0.6, linestyle="dashed")
    ax.plot(np.ones(num_lines) * x[i], y, zs=z[i], color='black', alpha=0.6, linestyle="dashed")

# Add arrows at random intersections
for i in range(num_lines):
    for j in range(num_lines):
        # Get a random direction
        direction = directions[i, j]
        # Get coordinates for the arrow based on the random direction
        arrow_coords = vertices_relative[direction]
        # Plot the arrow at the intersection point
        ax.quiver(x[i], y[j], z[i], arrow_coords[0], arrow_coords[1], arrow_coords[2], length=0.12, normalize=True)

# Set plot labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_box_aspect([1,1,0.175])
ax.set_axis_off()
ax.set_facecolor("white")

for cc in range(4):
    for dd in range(4):
        vector1 = vertices_relative[cc]
        vector2 = vertices_relative[dd]
        angle_degrees = np.degrees(np.arccos(np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))))
        print(angle_degrees)

plt.show()