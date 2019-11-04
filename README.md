# CIVL556_Hull
Term project for CIVL 556 W1 2019.
Written by Katherine Westerlund

Program intakes parameters for a given number of bounding boxes for subsystems within the submarine (in x,y coordinates
in meters).
Also intakes the maximum length of the hull, and the maximum width (both in meters).

Program takes these parameters and runs the "<<INSERT>>" optimization algorithm to minimize the sum of the volume of
the hull and the length of the hull.

Constraints on the solution space are:
> the Length Constraint (hull length must be less than the given maximum)
    > accomplished by comparing the hull length of the iteration to the given maximum
> the Width Constraint (hull width must be less than the given maximum)
    >  accomplished by comparing the hull width of the iteration to the given maximum
> the Bounding Box constraints (at every bounding point given, the point must be contained within the hull)
    > accomplished by calculating the height of the hull at the axial location of every bounding point, and comparing
    this height to the height of the bounding constraint.

Solution Reporting
> the solution is displayed, TODO: along with a plot of the bounding points and the final hullform
> TODO: sensitivity reporting on the optimal solution (how close are we to cutting off feets)

2D System mapping
1. intake units?
2. intake bounding dimensions and locations along the submarine body
3. intake other bounding considerations (eg. length)
4. intake airfoils for consideration (one at a time perform:)
    a. scale body to meet all bounding considerations
    b. report Cd
    c. calculate hull volume (assuming rotational integration)
    d. report values and store
5. loop through 4 until all data collected and reported
6. return optimal airfoil for human consideration, along with hull parameters
