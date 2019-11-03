# CIVIL 556 FINAL PROJECT
# Optimization Modelling and hull form optimization program
# Written by Katherine Westerlund
# Written 2019-09-25
# Updates written 2019-11-02

import numpy as np
from scipy.optimize import minimize


def optimization_function(length, thickness):
    # attempting to make the shortest submarine possible that can contain all the boxes
    optimum = volume_revolution(length, thickness) + length
    return optimum


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

    max_length = input("Please enter the maximum submarine length (in meters) ")

    max_width = input("Please enter the maximum width (in meters) ")

    return bounding_x, bounding_y, max_length, max_width, bounding_number


def height_constraint(constraints_x, constraints_y, bounding_number, thickness):
    # for every combo of x and y limits, check whether the constraint points are contained within the hull form
    index = 0
    while index > bounding_number:
        t = airfoil_at_point(thickness, constraints_x[index])
        if t < constraints_y[index]:
            return 0
        index = index + 1
    return 1


def length_constraint(max_length, length):
    if max_length >= length > 0:
        return 1
    else:
        return 0


def airfoil_at_point(t, x):
    y_t = 5*t*(0.2969*x**(1/2) - 0.126*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    return float(y_t)


def reporting(thickness, length):
    # make it print a pretty graph
    # also interested in calling the volume_revolution function

    return


def volume_revolution(thickness, length):
    # write a function to populate a data table with ~100 points along the airfoil
    count = 0
    data_points = []
    scale = length/1.0
    while count < length:
        data_points = data_points.insert(count, scale*airfoil_at_point(thickness, count))
        count = count + length/100
    volume = np.trapz(data_points)

    return volume


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
x0 = np.array(length, thickness)

# bounding total solution area
bound_horizontal = max_length
bound_vertical = max_width

bnds = np.array(bound_horizontal, bound_vertical)

# CONSTRAINTS
# con1 ensures bounding points are contained within the hull form
con1 = {'type': 'ineq', 'fun': height_constraint}

# con2 ensures the maximum length of the submarine is shorter than the maximum length
con2 = {'type': 'ineq', 'fun': length_constraint}

# array of constraints to be passed to the optimization function
cons = np.typeDict(con1, con2)

solution = minimize(optimization_function(length, thickness), x0, method='SLSQP', bounds=bnds, constraints=cons)
