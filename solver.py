from pyibex import Interval, Function, IntervalVector, CtcFwdBwd

from iodgsol import *

# Read the DGsol file and convert it into a list of constraints
constraints = read_dgsol_file('dgsol\data\data_set_1\graph.01.data')
print('Number of constraints: {}'.format(len(constraints)))

# Find the maximum point ID in the constraints
max_id = find_max_point_id(constraints)
print('Maximum point ID: {}'.format(max_id))


n = 3
border = 100


# Remove constraints with point IDs greater than n
constraints = remove_constraints_with_high_id(constraints, n)
print('Number of constraints after filtering: {}'.format(len(constraints)))

x, y, z = [Interval(-border,border) for _ in range(n )], [Interval(-border,border) for _ in range(n )], [Interval(-border,border) for _ in range(n )]


ctc = None

#distance constraints
for constraint in constraints:
    point1, point2 = constraint["points"]
    f = Function("x[{}]".format(n), "y[{}]".format(n), "z[{}]".format(n), "(x[{}]-x[{}])^2+(y[{}]-y[{}])^2+(z[{}]-z[{}])^2".format(point1-1,point2-1,point1-1,point2-1,point1-1,point2-1))
    ctc_f = CtcFwdBwd(f, Interval(constraint["dmin"]**2, constraint["dmax"]**2))
    if ctc is None:
        ctc = ctc_f
    else:
        ctc = ctc & ctc_f

# triangle inequality constraints
from itertools import combinations
for A, B, C in combinations(range(n), 3):
    # Constraint 1: AB + BC >= AC
    
    f1 = Function("x[{}]".format(n), "y[{}]".format(n), "z[{}]".format(n), 
                  "sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) + sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) - sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2)".format(A, B, A, B, A, B, B, C, B, C, B, C, A, C, A, C, A, C))
    ctc_f1 = CtcFwdBwd(f1, Interval(0, float('inf')))
    ctc = ctc & ctc_f1
    # Constraint 2: AB + AC >= BC
    f2 = Function("y[{}]".format(n), "x[{}]".format(n),  "z[{}]".format(n), 
                  "sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) + sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) - sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2)".format(A, B, A, B, A, B, B, C, B, C, B, C, A, C, A, C, A, C))
    ctc_f2 = CtcFwdBwd(f2, Interval(0, float('inf')))
    ctc = ctc & ctc_f2
    # Constraint 3: BC + AC >= AB
    f3 = Function("x[{}]".format(n),  "z[{}]".format(n), "y[{}]".format(n),
                  "sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) + sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2) - sqrt((x[{}]-x[{}])^2 + (y[{}]-y[{}])^2 + (z[{}]-z[{}])^2)".format(A, B, A, B, A, B, B, C, B, C, B, C, A, C, A, C, A, C))
    ctc_f3 = CtcFwdBwd(f3, Interval(0, float('inf')))
    ctc = ctc & ctc_f3


box = [Interval(-border,border) for _ in range(3*n )]
box[0] = Interval(0,0)
box[1] = Interval(0,0)
box[2] = Interval(0,0)
box[3] = Interval(0,border )
box[4] = Interval(0,0)
box[5] = Interval(0,0)
box[7] = Interval(0,border)
box[8] = Interval(0,0)
initial_box = IntervalVector(box)


stack = [initial_box]
stack_acc = []
stack_rej = []
stack_unc = []

tau = 0.02

while stack:
    box = stack.pop()
    original_box = box.copy()

    # Apply the contractor
    ctc.contract(box)


    # All points of the box are inside the constraints
    if (original_box == box) and (ctc.contract(box) != None):
        stack_acc.append(original_box)
    # No points of the box are inside the constraints
    elif box.is_empty():
        stack_rej.append(original_box)
    # Size of the box is big enough for further splitting
    elif box.max_diam() > tau:
        idx = max(range(len(box)), key=lambda i: box[i].diam())
        box1, box2 = box.bisect(idx)
        stack.append(box1)
        stack.append(box2)
    else:
        stack_unc.append(original_box)



print(len(stack_acc), len(stack_rej), len(stack_unc)) 

#print(stack_unc[0])



for su in stack_unc:
    cloud = []
    for i in range(n):
        cloud.append((su[3*i].mid(), su[3*i+1].mid(), su[3*i+2].mid()))
    #print(cloud, constraints)
    if check_constraints(cloud, constraints):
        print("found")
        print(cloud)
        break
    #print(check_constraints(cloud, constraints))


import matplotlib.pyplot as plt
import numpy as np

# Function to plot a circle
def plot_circle(center, radius,  label):
    theta = np.linspace(0, 2*np.pi, 100)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    plt.plot(x, y,  label=label)

# Function to plot a rectangle (box)
def plot_rectangle(lower_left, upper_right, color):
    x1, y1 = lower_left
    x2, y2 = upper_right
    plt.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], color=color)
# Function to plot a point
def plot_point(point,  label):
    x, y = point
    plt.scatter(x, y,  zorder=5, label=label)

# Plot the constraints
for i,(x,y,z) in enumerate(cloud):
    plot_point((x,y),  label='Point {}'.format(i+1))

for constraint in constraints:
    x,y,z = cloud[constraint["points"][0]-1]
    plot_circle((x,y), constraint["dmin"], label='dmin {}{}'.format(constraint["points"][0], constraint["points"][1]))
    plot_circle((x,y), constraint["dmax"], label='dmax {}{}'.format(constraint["points"][0], constraint["points"][1]))

for su in stack_unc:
    for i in range(1, n):
        x1, y1, z1 = su[3*i][0], su[3*i+1][0], su[3*i+2][0]
        x2, y2, z2 = su[3*i][1], su[3*i+1][1], su[3*i+2][1]
        plot_rectangle((x1, y1), (x2, y2), 'red')

# Additional plot settings
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.legend()
plt.xlabel("x")
plt.ylabel("y")

plt.show()