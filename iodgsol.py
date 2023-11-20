import numpy as np
def read_dgsol_file(file_path):
    """
    Reads a DGsol file and converts it into a list of distance constraints.
    
    Each line in the DGsol file is expected to represent a distance constraint between two points,
    including the minimum and maximum allowable distances. This function parses each line and 
    converts it into a dictionary format suitable for further processing or constraint handling.
    
    Args:
        file_path (str): The path to the DGsol file.
    
    Returns:
        list of dict: A list of distance constraints where each constraint is represented as a dictionary.
                      Each dictionary has the following keys:
                        - 'points': a list of two integers representing the points involved in the constraint.
                        - 'dmin': a float representing the minimum allowable distance between the points.
                        - 'dmax': a float representing the maximum allowable distance between the points.
                      
    Example of a returned constraint:
        {
            'points': [1, 2],
            'dmin': 0.5,
            'dmax': 1.5
        }
    """
    constraints = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.split()
            point1 = int(parts[0])
            point2 = int(parts[1])
            distance_min = float(parts[2].replace('E', 'e'))  # Convert to float
            distance_max = float(parts[3].replace('E', 'e'))  # Convert to float
            constraints.append({"points": [point1, point2], "dmin": distance_min, "dmax": distance_max})
                
    return constraints

def remove_constraints_with_high_id(constraints, n):
    """
    Removes constraints from the list where any point's ID is greater than n.

    Args:
        constraints (list of dict): The list of constraints where each constraint is a dictionary.
                                    Each dictionary must contain a key 'points' with a list of two point IDs.
        n (int): The maximum allowable point ID.

    Returns:
        list of dict: A filtered list of constraints where all points have IDs less than or equal to n.
    """
    filtered_constraints = [constraint for constraint in constraints 
                            if constraint["points"][0] <= n and constraint["points"][1] <= n]
    return filtered_constraints

def find_max_point_id(constraints):
    """
    Finds the maximum point ID from a list of constraints.

    Args:
        constraints (list of dict): A list of constraints where each constraint is represented as a dictionary.
                                    Each dictionary should contain a key 'points' with a list of two point IDs.

    Returns:
        int: The maximum point ID found in the constraints.
    """
    max_id = 0
    for constraint in constraints:
        max_id = max(max_id, max(constraint["points"]))
    return max_id

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw_3d_point_cloud(points):
    """
    Draws a 3D cloud of points with labels.

    Args:
        points (list of tuples): A list of points, where each point is a tuple (x, y, z).
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for i, (x, y, z) in enumerate(points):
        ax.scatter(x, y, z, marker='o')
        ax.text(x, y, z, f'{i}', color='blue')

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()

import math

def check_constraints(points, constraints):
    """
    Checks if a cloud of points satisfies a list of constraints.

    Args:
        points (dict): A dictionary of points where the key is the point ID and the value is a tuple (x, y, z) representing the point's coordinates.
        constraints (list of dict): A list of constraints where each constraint is represented as a dictionary.
                                    Each dictionary should contain a key 'points' with a list of two point IDs, 'dmin' and 'dmax' for minimum and maximum distance.

    Returns:
        bool: True if all constraints are satisfied, False otherwise.
    """
    def euclidean_distance(point1, point2):
        return math.sqrt(sum((p1 - p2) ** 2 for p1, p2 in zip(point1, point2)))

    for constraint in constraints:
        point1, point2 = constraint["points"]
        actual_distance = euclidean_distance(points[point1-1], points[point2-1])
        if not (constraint["dmin"] <= actual_distance <= constraint["dmax"]):
            return False

    return True

def draw_constraints_as_spheres(points, constraints):
    """
    Draws a 3D plot with points and spheres representing the constraints.

    Args:
        points (dict): A dictionary of points with point ID as keys and (x, y, z) tuples as values.
        constraints (list of dict): A list of constraints, each with 'points', 'dmin', and 'dmax' keys.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Draw points
    for point_id, (x, y, z) in enumerate(points):
        ax.scatter(x, y, z, color='blue')
        ax.text(x, y, z, f'{point_id}')

    # Function to draw a sphere
    def draw_sphere(x_center, y_center, z_center, radius, color):
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x_sphere = radius * np.cos(u) * np.sin(v) + x_center
        y_sphere = radius * np.sin(u) * np.sin(v) + y_center
        z_sphere = radius * np.cos(v) + z_center
        ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color=color, alpha=0.3)

    # Draw spheres for constraints
    for constraint in constraints:
        point1_id, point2_id = constraint['points']
        x1, y1, z1 = points[point1_id]


        # Draw spheres for the minimum and maximum distances
        draw_sphere(x1, y1, z1, constraint['dmin'], 'green')
        draw_sphere(x1, y1, z1, constraint['dmax'], 'red')


    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()




if __name__ == '__main__': 
    # Read the DGsol file and convert it into a list of constraints
    constraints = read_dgsol_file('dgsol\data\data_set_1\graph.01.data')
    print('Number of constraints: {}'.format(len(constraints)))

    # Find the maximum point ID in the constraints
    max_id = find_max_point_id(constraints)
    print('Maximum point ID: {}'.format(max_id))

    # Remove constraints with point IDs greater than 5
    constraints = remove_constraints_with_high_id(constraints, 5)
    print('Number of constraints after filtering: {}'.format(len(constraints)))
    print('Constraints: {}'.format(constraints))

    # # Example usage
    # points = {1: (1, 2, 3), 2: (4, 5, 6), 3: (7, 8, 9)}  # Replace with your points
    # constraints = [{"points": [1, 2], "dmin": 2, "dmax": 5},
    #             {"points": [2, 3], "dmin": 3, "dmax": 7}]  # Replace with your constraints
    # result = check_constraints(points, constraints)
    # print('Do the points satisfy the constraints?', result)


    # # Example usage
    # points = [(1, 2, 3), (4, 5, 6), (7, 8, 9)]  # Replace with your list of points
    # draw_3d_point_cloud(points)

    # Example usage
    points = {1: (0, 0, 0), 2: (1, 1, 1), 3: (2, 2, 2)}
    constraints = [{"points": [1, 2], "dmin": 0.5, "dmax": 1.5},
                {"points": [2, 3], "dmin": 1, "dmax": 2}]
    draw_constraints_as_spheres(points, constraints)