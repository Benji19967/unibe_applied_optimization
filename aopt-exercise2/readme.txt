Theoretical Exercises (1.0 pt)

The proofs are based on methods extracted from course material (slides) and specific methods are detailed in the solutions.

Done by: Horacio Lisdero

Programming Exercise: Convexity Test (1 pt)

1. Solution Description

The algorithm proceeds as follows:
- Outer loop samples 1,000,000 random pairs of points (x, y).
- Inner loop check (n_evals) uniformly spaced points on the segment between x and y.
- If a counterexample is found, the function immediately prints the error details and returns false.
- If all trials pass, the function returns true (likely convex).

2. Encountered Problems and Solutions

- Numerical Precision: Due to floating-point errors, standard comparison had a risk of producing false negatives.
  - Solution: Introduced epsilon = 1e-9: LHS > RHS + epsilon.
- Theta Sampling: Ensuring uniform sampling including endpoints (t=0 and t=1) across n_evals steps.
  - Solution: Calculated t as k / (n_evals - 1.0) for k = 0 to n_evals - 1.

3. Done by: Beñat Solaegui
