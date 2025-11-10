# Exercise 6 - Newton's Method

## Theory: Affine Invariance (1 + Bonus 0.5 pt)

Done by: Beñat Solaegui

Solved by applying the ideas from the lecture slides.

## Programming Exercise (2pt)

### Standard Newton Method

Done by: Benjamin Labrecque

Implement Newton's method as described in the algorithm from the lecture slides. 

### Newton method with modified hessian

Done by: Benjamin Labrecque

Same as previous algorithm but update the Hessian iteratively if the Cholesky 
decomposition fails.

### Newton Methods vs Gradient descent

Done by: Horacio Lisdero

For a timing report in table format, see `assets/summary.csv`.

From `summary.csv` we see that Newton's method in any of its variants reaches the optimum significantly faster than gradient descent.

A more detailed timing report is in `assets/TIMING.md`. We can see that gradient descent performs more evaluations of the function.

The screenshots can be found in `assets/screenshots`.

We can see that the optimization results look the same, regardless of the method used,
except for `spring_noLength_gradient_descent_20x20.png` where gradient descent failed.

### Bonus (0.5 pt)

Done by: Horacio Lisdero

Implement the projection following the idea in "Strategy I: Direct Eigenvalue Modification" in the lecture slides.

Executed the standard Newton method and added it to the table.

We observe in the report that:

- Combining the spring element with the local positive hessian and the Newton's method with projected Hessian
gives similar results as using the standard method with the "fix" in the spring element.
The reason is that the choice of the method overrides the "fix" in the element.
