# Exercise 7 - L-BFGS and Gauss Newton Method

## Theory: Secant Equation (0.5 pt) + Curvature Condition (0.5 pt)

Done by: Horacio Lisdero

Solved by applying the ideas from the lecture slides.

## Programming (2 pts)

### L-BFGS

Done by:

### Gauss-Newton Method

Done by:

### Comparison of different methods

Done by: Benjamin Labrecque

Created bash script (compare_methods.sh) to generate timings. 
Created python script (compare_timings.py) to parse num iterations, total run times, 
and objective function value progression from timing outputs. Also created plots for the
metrics -- see:
 - assets/num_iterations
 - assets/obj_function_values
 - assets/total_run_times

Fastest: projected newton seems to consistently perform best.

A summary of the metrics can be found at assets/summary.csv

### A line search algorithm for the Wolfe conditions (Bonus (1 pt))

Done by:
