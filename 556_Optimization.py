# CIVIL 556 FINAL PROJECT
# Optimization Modelling and hull form optimization program
# Written by Katherine Westerlund
# Written 2019-09-25
# Updates written 2019-11-02

import numpy as np
from scipy.optimize import minimize


def intake_data():

    # To Do: Display model for context for bounding data intake
    # intake number of bounding boxes
    # self.is_not_used()

    bounding_number = input('How many bounding boxes do you need? ')
    bounding_number = float(bounding_number)

    if isinstance(bounding_number, float) == 0:
        return print("Error, number of bounding boxes is not numerical")

    bounding_x = [0]
    bounding_y = [0]

    count = 0

    while count < bounding_number:
        x = float(input("Please enter x-component (in meters) of data point " + str(count)))
        bounding_x.insert(count, x)

        y = float(input("Please enter the y-component (in meters) of data point " + str(count)))
        bounding_y.insert(count, y)

        count = count + 1

    count = 0
    print("The bounding points (in meters) collected are: \n")
    while count < bounding_number:
        print("(" + str(bounding_x[count]) + " " + str(bounding_y[count]) + ")")
        count = count + 1

    max_l = input("Please enter the maximum submarine length (in meters) ")

    max_w = input("Please enter the maximum width (in meters) ")

    return bounding_x, bounding_y, max_l, max_w, bounding_number


def height_constraint(x):
    # for every combo of x and y limits, check whether the constraint points are contained within the hull form
    # front_gap modifies x distance between the airfoil front and the x-component of the first bounding pt
    cons_x = constraints_x
    cons_y = constraints_y
    length_testing = x[0]
    front_gap = x[2]
    index = 0
    t = x[1]
    while index < bounding_number:
        height = airfoil_at_point(t, (cons_x[index] + front_gap), length_testing)
        if height < cons_y[index]:
            return -1
        index = index + 1
    return 1


def length_constraint(x):
    max_l = float(max_length)
    length = float(x[0])
    if max_l >= length:
        return 1
    else:
        return 0


def airfoil_at_point(t, point, length):
    # input real length, real t, real x
    scale = length/1.0
    # convert to x = [0,1.0] by dividing by the real length
    point = point/length
    # scale down thickness to the theoretical thickness of a 1m sub by dividing by the scale
    t = t/scale
    # t is the max thickness, x is the position from 0-1.0 along the airfoil
    y_t = 5*t*(0.2969*point**(1/2) - 0.126*point - 0.3516*point**2 + 0.2843*point**3 - 0.1015*point**4)
    # convert based on scale of length
    y_t = y_t*scale
    return float(y_t)


def volume_revolution(x):
    # write a function to populate a data table with ~200 points along the airfoil
    test_len = x[0]
    hold = 0
    data_points = np.zeros(shape=(1, 1))
    x_components = np.zeros(shape=(1, 1))
    while hold <= test_len:
        data_points = np.append(data_points, (airfoil_at_point(x[1], hold, x[0])))
        x_components = np.append(x_components, hold)
        hold = hold + test_len/100
    # take the trapezoidal integral of ~200 data points evenly spaced along the real length of the submarine
    print(data_points)
    integral = np.trapz(data_points, dx=test_len/100, axis=0)
    print(integral)
    # calculate the rotational volume from the integral
    volume = integral*np.pi
    print(volume)
    return volume


def optimization_function(x):
    # attempting to make the shortest submarine possible that can contain all the boxes
    # implement coefficients in the future
    print("Testing optimization function")
    print(x)
    optimum = x[0]
    return optimum


def reporting(x):
    # make it print a pretty graph
    # also interested in calling the volume_revolution function
    # display 2D curve of optimum hull form with bounding points placed on it
    print(x)
    return


print("Welcome to the CIVIL 556 Hull form optimization modeller! \n "
      "Please follow the instructions for data entry, and enjoy!")

# MODEL PARAMETERS,
constraints_x, constraints_y, max_length, max_width, bounding_number = intake_data()

# INITIALIZING SOLUTION VARIABLES
# initialized length
x0_length = 2

# thickness at max point
x0_thickness = 0.2

# front_space is the difference between the start of the airfoil and the first x component of bounding
x0_front_space = 0.1

# initial solution array
x0 = np.array([x0_length, x0_thickness, x0_front_space])

# BOUNDING
# bounding total solution area
bound_horizontal = [1.2, max_length]
bound_vertical = [0.1, max_width]
bound_front = (0.01, 1)

# bounding array for each variable
bnds = np.array([bound_horizontal, bound_vertical, bound_front])

# CONSTRAINTS
# con1 ensures bounding points are contained within the hull form
con1 = {'type': 'ineq', 'fun': height_constraint}

# con2 ensures the maximum length of the submarine is shorter than the maximum length
con2 = {'type': 'ineq', 'fun': length_constraint}

# array of constraints to be passed to the optimization function
cons = con1, con2

solution = minimize(optimization_function, x0, method='SLSQP', bounds=bnds, constraints=con1)

x = solution.x
print("Solution")
print(solution)
print(x[0])
print(x[1])
print(x[2])
