# CIVL556_Hull
Term project for CIVL 556 2019. Shape optimization algorithm for optimizing shape of submarine hullforms.

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
    a. would also be good to report sensitivity of design to changes in bounding considerations (how close are we to cutting off feets)

