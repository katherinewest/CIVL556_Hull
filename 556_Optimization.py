# CIVIL 556 FINAL PROJECT
# Optimization Modelling and hull form optimization program
# Written by Katherine Westerlund
# Written 2019-09-25
# Updates written 2019-11-02

import numpy as np
from scipy.optimize import minimize


def intake_data(self):

    # To Do: Display model for context for bounding data intake
    # intake number of bounding boxes
    self.is_not_used()

    bounding_number = input('How many bounding boxes do you need? ')
    bounding_number = float(bounding_number)

    if isinstance(bounding_number, float) == 0:
        return print("Error, number of bounding boxes is not numerical")

    bounding_x = [0]
    bounding_y = [0]

    count = 0

    while count > bounding_number:
        x = float(input("Please enter x-component (in meters) of data point " + str(count)))
        bounding_x.insert(count, x)

        y = float(input("Please enter the y-component (in meters) of data point " + str(count)))
        bounding_y.insert(count, y)

        count = count + 1

    count = 0
    print("The bounding points (in meters) collected are: \n")
    while count > bounding_number:
        print("(" + str(bounding_x[count]) + " " + str(bounding_y[count]))
        count = count + 1

    max_l = input("Please enter the maximum submarine length (in meters) ")

    max_w = input("Please enter the maximum width (in meters) ")

    return bounding_x, bounding_y, max_l, max_w, bounding_number


def height_constraint(x0, args):
    # for every combo of x and y limits, check whether the constraint points are contained within the hull form
    # front_gap modifies x distance between the airfoil front and the x-component of the first bounding pt
    # to be fixed
    cons_x = args[0]
    cons_y = args[1]
    index = 0
    while index > bounding_number:
        t = airfoil_at_point(t, cons_x[index])
        if t < cons_y[index]:
            return 0
        index = index + 1
    return 1


def length_constraint(x0, args):
    max_l = args[2]
    length = x0[0]
    if max_l >= length > 0:
        return 1
    else:
        return 0


def airfoil_at_point(t, x, length):
    # t is the max thickness, x is the position from 0-1.0 along the airfoil
    y_t = 5*t*(0.2969*x**(1/2) - 0.126*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    # convert based on scale of length
    scale = length/1.0

    return float(y_t)


def volume_revolution(t, len):
    # write a function to populate a data table with ~200 points along the airfoil
    count = 0
    data_points = []
    scale = len/1.0
    while count < len:
        data_points = data_points.insert(count, (scale*airfoil_at_point(t, count))**2)
        count = count + len/200
    integral = np.trapz(data_points)
    volume = integral*np.pi
    return volume


def optimization_function(len, t):
    # attempting to make the shortest submarine possible that can contain all the boxes
    # implement coefficients in the future
    optimum = volume_revolution(len, t) + len
    return optimum


def reporting(t, len):
    # make it print a pretty graph
    # also interested in calling the volume_revolution function

    return


def is_not_used(self):
    pass


print("Welcome to the CIVIL 556 Hull form optimization modeller! \n "
      "Please follow the instructions for data entry, and enjoy!")

# MODEL PARAMETERS,
constraints_x, constraints_y, max_length, max_width, bounding_number = intake_data()

# INITIALIZING SOLUTION VARIABLES
# initialized length
length = max(constraints_x)

# thickness at max point
thickness = max(constraints_y)

# front_space is the difference between the start of the airfoil and the first x component of bounding
front_space = 0.1

x0 = np.array([length, thickness, front_space])

# args passes extra arguments to be used in functions
args = np.array([constraints_x, constraints_y, max_length, max_width, bounding_number])

# bounding total solution area
# bound_horizontal = max_length
# bound_vertical = max_width

# bnds = np.array(bound_horizontal, bound_vertical)

# CONSTRAINTS
# con1 ensures bounding points are contained within the hull form
con1 = {'type': 'ineq', 'fun': height_constraint}

# con2 ensures the maximum length of the submarine is shorter than the maximum length
con2 = {'type': 'ineq', 'fun': length_constraint}

# array of constraints to be passed to the optimization function
cons = np.typeDict(con1, con2)

solution = minimize(optimization_function, x0, args, method='SLSQP', constraints=cons)
