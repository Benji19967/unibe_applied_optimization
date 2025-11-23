# Applied Optimization - University of Bern

Exercises for course applied optimization at the University of Bern

## Getting started with an exercise

```shell
cd aopt-exercise<n>
mkdir build
cd build
CMAKE_POLICY_VERSION_MINIMUM=3.5 cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make
```

## Selecting a subset of tests to run

```shell
cd build
ctest -R <regex to select tests>
```

## Visualization tools

### Contour and Surface plots

https://slides.cgg.unibe.ch/aopt20/plots/plot2d_create.html

Using a custom csv file:
https://slides.cgg.unibe.ch/aopt20/plots/plot2d_custom.html

### Spring System

https://slides.cgg.unibe.ch/aopt20/plots/springsys.html

### Taylor approx and descent methods

https://cgg.unibe.ch/interactive-education/

- Taylor polynomial approx in 1-D and 2-D
- Descent methods in 2-D

Credits to Martin Heistermann
