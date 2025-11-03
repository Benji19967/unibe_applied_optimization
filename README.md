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

## Interactive visualization tools

https://cgg.unibe.ch/interactive-education/

Credits to Martin Heistermann
