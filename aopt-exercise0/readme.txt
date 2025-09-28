# Exercise 0

## Implementing objective functions (1 pt)

Done by: Benjamin Labrecque

Solved by direct implementation of the function definitions.

### Screenshots

![3D surface of 2D function](./assets/func2d_0.png "3D surface of 2D function")

![Isocontours of the non-convex function](./assets/funcNonConvex_0.png "Isocontours of the non-convex function")

## Grid Search (1 pt)

Done by: Horacio Lisdero

In the 2-D case, solved by using two nested for loops.

In the N-D case, solved by an iterative approach that treats the n-dimensional index as
a "number" with N digits, where each digit is incremented and then carried over to the next.

All the N-D points at which the function is to be evaluated can be expressed by the equation:

```c
x = lower + index * delta;
```

Where `lower` is a N-D vector with the coordinates of the lower corner of the cuboid, `index` is the current N-D index and `delta` is an N-D vector with the step size to take for every dimension. Note: `*` represents the Hadamard product (coefficient-wise product) in this case.

### How it was solved?

## Theory Background Exercises (optional)

Done by: Beñat Solaegui
