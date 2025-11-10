# Timing report

For an easier to view format, see `summary.csv`.

Methods:

- Standard Newton's method
- Newton's method with projected Hessian
- Gradient descent

## Standard Newton's method

### No length

```txt
Initial MassSpring system energy is 2.26635e+07
Saving initial spring graph to spring_noLength_standard_5x5_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.000261s
total time evaluation : 9.2e-05s  (35.249 %)
eval_f time   : 0.00000s  ( #evals: 4 -> avg 0.00000s )
eval_grad time: 0.00000s  ( #evals: 2 -> avg 0.00000s, factor: 6.00000)
eval_hess time: 0.00009s  ( #evals: 2 -> avg 0.00004s, factor: 176.00000)
Saving optimized spring graph to spring_noLength_standard_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.30682e+07
Saving initial spring graph to spring_noLength_standard_10x10_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.001125s
total time evaluation : 0.000378s  (33.6 %)
eval_f time   : 0.00002s  ( #evals: 7 -> avg 0.00000s )
eval_grad time: 0.00002s  ( #evals: 3 -> avg 0.00001s, factor: 2.18750)
eval_hess time: 0.00035s  ( #evals: 3 -> avg 0.00012s, factor: 50.60417)
Saving optimized spring graph to spring_noLength_standard_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.05053e+08
Saving initial spring graph to spring_noLength_standard_20x20_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.005792s
total time evaluation : 0.001616s  (27.9006 %)
eval_f time   : 0.00008s  ( #evals: 12 -> avg 0.00001s )
eval_grad time: 0.00007s  ( #evals: 4 -> avg 0.00002s, factor: 2.55556)
eval_hess time: 0.00147s  ( #evals: 4 -> avg 0.00037s, factor: 54.29630)
Saving optimized spring graph to spring_noLength_standard_20x20_opt_*.csv
```

### With length

```txt
Initial MassSpring system energy is 2.43845e+07
Saving initial spring graph to spring_withLength_standard_5x5_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.000605s
total time evaluation : 0.00018s  (29.7521 %)
eval_f time   : 0.00000s  ( #evals: 20 -> avg 0.00000s )
eval_grad time: 0.00001s  ( #evals: 10 -> avg 0.00000s, factor: 10.00000)
eval_hess time: 0.00017s  ( #evals: 10 -> avg 0.00002s, factor: 348.00000)
Saving optimized spring graph to spring_withLength_standard_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.90827e+07
Saving initial spring graph to spring_withLength_standard_10x10_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.003002s
total time evaluation : 0.000969s  (32.2785 %)
eval_f time   : 0.00004s  ( #evals: 22 -> avg 0.00000s )
eval_grad time: 0.00004s  ( #evals: 11 -> avg 0.00000s, factor: 1.63636)
eval_hess time: 0.00089s  ( #evals: 11 -> avg 0.00008s, factor: 40.40909)
Saving optimized spring graph to spring_withLength_standard_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.29333e+08
Saving initial spring graph to spring_withLength_standard_20x20_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.017239s
total time evaluation : 0.003942s  (22.8668 %)
eval_f time   : 0.00022s  ( #evals: 26 -> avg 0.00001s )
eval_grad time: 0.00017s  ( #evals: 12 -> avg 0.00001s, factor: 1.69820)
eval_hess time: 0.00355s  ( #evals: 12 -> avg 0.00030s, factor: 34.60811)
Saving optimized spring graph to spring_withLength_standard_20x20_opt_*.csv
```

### With length with positive local Hessian

```txt
Initial MassSpring system energy is 2.43845e+07
Saving initial spring graph to spring_withLengthPositiveLocalHessian_standard_5x5_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.00181s
total time evaluation : 0.001285s  (70.9945 %)
eval_f time   : 0.00000s  ( #evals: 23 -> avg 0.00000s )
eval_grad time: 0.00001s  ( #evals: 11 -> avg 0.00000s, factor: 23.00000)
eval_hess time: 0.00127s  ( #evals: 11 -> avg 0.00012s, factor: 2661.72727)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_standard_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.90827e+07
Saving initial spring graph to spring_withLengthPositiveLocalHessian_standard_10x10_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.00735s
total time evaluation : 0.005133s  (69.8367 %)
eval_f time   : 0.00005s  ( #evals: 24 -> avg 0.00000s )
eval_grad time: 0.00004s  ( #evals: 12 -> avg 0.00000s, factor: 1.42308)
eval_hess time: 0.00504s  ( #evals: 12 -> avg 0.00042s, factor: 194.00000)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_standard_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.29333e+08
Saving initial spring graph to spring_withLengthPositiveLocalHessian_standard_20x20_*.csv
******** Newton Method ********
######## Timing statistics ########
total time    : 0.038912s
total time evaluation : 0.023425s  (60.1999 %)
eval_f time   : 0.00026s  ( #evals: 26 -> avg 0.00001s )
eval_grad time: 0.00017s  ( #evals: 13 -> avg 0.00001s, factor: 1.34630)
eval_hess time: 0.02300s  ( #evals: 13 -> avg 0.00177s, factor: 178.94942)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_standard_20x20_opt_*.csv
```

## Newton's method with projected Hessian

### No length

```txt
Initial MassSpring system energy is 2.26635e+07
Saving initial spring graph to spring_noLength_projected_hessian_5x5_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.000196s
total time evaluation : 6.8e-05s  (34.6939 %)
eval_f time   : 0.00000s  ( #evals: 4 -> avg 0.00000s )
eval_grad time: 0.00000s  ( #evals: 2 -> avg 0.00000s, factor: inf)
eval_hess time: 0.00006s  ( #evals: 2 -> avg 0.00003s, factor: inf)
Saving optimized spring graph to spring_noLength_projected_hessian_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.30682e+07
Saving initial spring graph to spring_noLength_projected_hessian_10x10_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.000945s
total time evaluation : 0.000326s  (34.4974 %)
eval_f time   : 0.00001s  ( #evals: 7 -> avg 0.00000s )
eval_grad time: 0.00001s  ( #evals: 3 -> avg 0.00000s, factor: 4.08333)
eval_hess time: 0.00030s  ( #evals: 3 -> avg 0.00010s, factor: 88.66667)
Saving optimized spring graph to spring_noLength_projected_hessian_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.05053e+08
Saving initial spring graph to spring_noLength_projected_hessian_20x20_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.005927s
total time evaluation : 0.001631s  (27.5181 %)
eval_f time   : 0.00008s  ( #evals: 12 -> avg 0.00001s )
eval_grad time: 0.00007s  ( #evals: 4 -> avg 0.00002s, factor: 2.59756)
eval_hess time: 0.00148s  ( #evals: 4 -> avg 0.00037s, factor: 54.07317)
Saving optimized spring graph to spring_noLength_projected_hessian_20x20_opt_*.csv
```

### With length

```txt
Initial MassSpring system energy is 2.43845e+07
Saving initial spring graph to spring_withLength_projected_hessian_5x5_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.000691s
total time evaluation : 0.000226s  (32.7062 %)
eval_f time   : 0.00000s  ( #evals: 20 -> avg 0.00000s )
eval_grad time: 0.00001s  ( #evals: 10 -> avg 0.00000s, factor: 12.00000)
eval_hess time: 0.00022s  ( #evals: 10 -> avg 0.00002s, factor: 438.00000)
Saving optimized spring graph to spring_withLength_projected_hessian_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.90827e+07
Saving initial spring graph to spring_withLength_projected_hessian_10x10_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.003235s
total time evaluation : 0.001054s  (32.5811 %)
eval_f time   : 0.00004s  ( #evals: 22 -> avg 0.00000s )
eval_grad time: 0.00004s  ( #evals: 11 -> avg 0.00000s, factor: 1.63636)
eval_hess time: 0.00097s  ( #evals: 11 -> avg 0.00009s, factor: 44.27273)
Saving optimized spring graph to spring_withLength_projected_hessian_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.29333e+08
Saving initial spring graph to spring_withLength_projected_hessian_20x20_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.02281s
total time evaluation : 0.005515s  (24.178 %)
eval_f time   : 0.00032s  ( #evals: 26 -> avg 0.00001s )
eval_grad time: 0.00027s  ( #evals: 12 -> avg 0.00002s, factor: 1.87090)
eval_hess time: 0.00493s  ( #evals: 12 -> avg 0.00041s, factor: 33.89630)
Saving optimized spring graph to spring_withLength_projected_hessian_20x20_opt_*.csv
```

### With length with positive local Hessian

```txt
Initial MassSpring system energy is 2.43845e+07
Saving initial spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_5x5_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.001881s
total time evaluation : 0.001379s  (73.3121 %)
eval_f time   : 0.00000s  ( #evals: 23 -> avg 0.00000s )
eval_grad time: 0.00003s  ( #evals: 11 -> avg 0.00000s, factor: 52.27273)
eval_hess time: 0.00135s  ( #evals: 11 -> avg 0.00012s, factor: 2829.00000)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_5x5_opt_*.csv
```

```txt
Initial MassSpring system energy is 7.90827e+07
Saving initial spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_10x10_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.007376s
total time evaluation : 0.005138s  (69.6584 %)
eval_f time   : 0.00005s  ( #evals: 24 -> avg 0.00000s )
eval_grad time: 0.00004s  ( #evals: 12 -> avg 0.00000s, factor: 1.45098)
eval_hess time: 0.00505s  ( #evals: 12 -> avg 0.00042s, factor: 198.03922)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_10x10_opt_*.csv
```

```txt
Initial MassSpring system energy is 4.29333e+08
Saving initial spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_20x20_*.csv
******** Newton Method with projected hessian ********
######## Timing statistics ########
total time    : 0.037917s
total time evaluation : 0.023097s  (60.9146 %)
eval_f time   : 0.00022s  ( #evals: 26 -> avg 0.00001s )
eval_grad time: 0.00017s  ( #evals: 13 -> avg 0.00001s, factor: 1.54667)
eval_hess time: 0.02270s  ( #evals: 13 -> avg 0.00175s, factor: 201.76000)
Saving optimized spring graph to spring_withLengthPositiveLocalHessian_projected_hessian_20x20_opt_*.csv
```

## Gradient descent

### No length

```txt
######## Timing statistics ########
total time    : 0.288159s
total time evaluation : 0.014652s  (5.08469 %)
eval_f time   : 0.00428s  ( #evals: 357005 -> avg 0.00000s )
eval_grad time: 0.01037s  ( #evals: 10000 -> avg 0.00000s, factor: 86.51532)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_noLength_gradient_descent_5x51_opt_*.csv
```

```txt
######## Timing statistics ########
total time    : 0.909113s
total time evaluation : 0.537953s  (59.1734 %)
eval_f time   : 0.49422s  ( #evals: 356683 -> avg 0.00000s )
eval_grad time: 0.04374s  ( #evals: 10000 -> avg 0.00000s, factor: 3.15656)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_noLength_gradient_descent_10x101_opt_*.csv
```

```txt
######## Timing statistics ########
total time    : 3.23327s
total time evaluation : 2.75016s  (85.0581 %)
eval_f time   : 2.56924s  ( #evals: 355992 -> avg 0.00001s )
eval_grad time: 0.18092s  ( #evals: 10000 -> avg 0.00002s, factor: 2.50684)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_noLength_gradient_descent_20x201_opt_*.csv
```

### With length

```txt
######## Timing statistics ########
total time    : 0.029593s
total time evaluation : 0.001243s  (4.20032 %)
eval_f time   : 0.00074s  ( #evals: 36980 -> avg 0.00000s )
eval_grad time: 0.00050s  ( #evals: 1021 -> avg 0.00000s, factor: 24.61940)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_withLength_gradient_descent_5x51_opt_*.csv
```

```txt
######## Timing statistics ########
total time    : 0.405147s
total time evaluation : 0.32703s  (80.7189 %)
eval_f time   : 0.31333s  ( #evals: 132232 -> avg 0.00000s )
eval_grad time: 0.01370s  ( #evals: 3647 -> avg 0.00000s, factor: 1.58545)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_withLength_gradient_descent_10x101_opt_*.csv
```

```txt
######## Timing statistics ########
total time    : 3.73879s
total time evaluation : 3.34401s  (89.4411 %)
eval_f time   : 3.20044s  ( #evals: 363866 -> avg 0.00001s )
eval_grad time: 0.14357s  ( #evals: 10000 -> avg 0.00001s, factor: 1.63231)
eval_hess time: 0.00000s  ( #evals: 0 -> avg -nans, factor: -nan)
Saving optimized spring graph to spring_withLength_gradient_descent_20x201_opt_*.csv
```
