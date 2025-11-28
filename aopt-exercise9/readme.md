# Exercise 9 - Equality Constrained Optimization II

## Theory

Done by: 

Solved by applying the ideas from the lecture slides.

## Programming

### Infeasible Start Newton Method (1 pt)

Done by: Beñat Solaegui

Solved by implementing the algorithm described in the lecture slides. 

### Hybrid Newton Method (Bonus 0.5 pt)

Done by: Beñat Solaegui

The algorithm was not explicitly given in the slides. I implemented it like this: when the iterate is infeasible (Ax-b larger than epsilon), it performs an infeasible-step Newton update. Once the iterate is sufficiently feasible, it switches to a feasible-step Newton update.


Ran the executable `EqualityConstrainedNewtonInfeasibleStart` for the three methods (fesible start, infeasible start and hybrid).
See the images in `assets/screenshots` and time statistics in `assets/statistics.txt`.
