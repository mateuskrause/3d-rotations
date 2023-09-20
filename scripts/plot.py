import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sympy import symbols, cos, sin, Matrix

def cube_with_rotation(rotation_matrix):
    # Create a smaller figure with three 3D axes side by side
    fig = plt.figure(figsize=(8, 4))  # Adjust the figsize as needed
    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132, projection='3d')
    
    ax3 = fig.add_subplot(133, projection='3d')

    # Define the scale factor for the cube
    scale = 2.0

    # Define the coordinates of the vertices of the scaled cube centered at (2, 2, 2)
    center = [2, 2, 2]
    vertices = [
        [-scale/2 + center[0], -scale/2 + center[1], -scale/2 + center[2]],
        [scale/2 + center[0], -scale/2 + center[1], -scale/2 + center[2]],
        [scale/2 + center[0], scale/2 + center[1], -scale/2 + center[2]],
        [-scale/2 + center[0], scale/2 + center[1], -scale/2 + center[2]],
        [-scale/2 + center[0], -scale/2 + center[1], scale/2 + center[2]],
        [scale/2 + center[0], -scale/2 + center[1], scale/2 + center[2]],
        [scale/2 + center[0], scale/2 + center[1], scale/2 + center[2]],
        [-scale/2 + center[0], scale/2 + center[1], scale/2 + center[2]]
    ]

    # Define the edges of the scaled cube by connecting the vertices
    edges = [
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7]
    ]

    # Set the limits for the plot regions (all 4x4)
    ax1.set_xlim([0, 4])
    ax1.set_ylim([0, 4])
    ax1.set_zlim([0, 4])

    ax2.set_xlim([0, 4])
    ax2.set_ylim([0, 4])
    ax2.set_zlim([0, 4])

    ax3.set_xlim([0, 4])
    ax3.set_ylim([0, 4])
    ax3.set_zlim([0, 4])

    # Set equal aspect ratio for all axes
    ax1.set_box_aspect([1, 1, 1])
    ax2.set_box_aspect([1, 1, 1])
    ax3.set_box_aspect([1, 1, 1])

    # Extract the x, y, and z coordinates of the vertices
    x, y, z = zip(*vertices)

    # Plot the cube edges in the first subplot (original cube)
    for edge in edges:
        x_vals = [x[edge[0]], x[edge[1]]]
        y_vals = [y[edge[0]], y[edge[1]]]
        z_vals = [z[edge[0]], z[edge[1]]]
        ax1.plot(x_vals, y_vals, z_vals, color='b')

    # Apply the rotation matrix to the cube's vertices
    rotated_vertices = [(rotation_matrix * (Matrix(vertex) - Matrix(center)) + Matrix(center)).tolist() for vertex in vertices]

    # Extract the rotated x, y, and z coordinates of the vertices
    x_rotated, y_rotated, z_rotated = zip(*rotated_vertices)

    # Plot the cube edges in the second subplot (rotated cube)
    for edge in edges:
        x_vals = [x_rotated[edge[0]], x_rotated[edge[1]]]
        y_vals = [y_rotated[edge[0]], y_rotated[edge[1]]]
        z_vals = [z_rotated[edge[0]], z_rotated[edge[1]]]
        ax2.plot(x_vals, y_vals, z_vals, color='r')

    # Plot both the original and rotated cubes in the third subplot
    for edge in edges:
        x_vals = [x[edge[0]], x[edge[1]]]
        y_vals = [y[edge[0]], y[edge[1]]]
        z_vals = [z[edge[0]], z[edge[1]]]
        ax3.plot(x_vals, y_vals, z_vals, color='b')

    for edge in edges:
        x_vals = [x_rotated[edge[0]], x_rotated[edge[1]]]
        y_vals = [y_rotated[edge[0]], y_rotated[edge[1]]]
        z_vals = [z_rotated[edge[0]], z_rotated[edge[1]]]
        ax3.plot(x_vals, y_vals, z_vals, color='r')

    # Show the plot
    plt.show()