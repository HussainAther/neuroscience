# Recreating weak pairwise correlations

This repository contains the code and output in recreating the work of Schneidman et al., (2006). 'Weak pairwise correlations imply strongly correlated network states in a neural population.' _Nature_.

The scripts use code by the [maxent_toolbox](https://github.com/orimaoz/maxent_toolbox).

* Data in `data`.
* Code in `src`.
* Output in `doc`.

## Working with small distributions of neurons (exhaustively)
Used in `smalldist.m`.
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/smalldist.png)

## Working with larger distributions of neurons (MCMC)
`largedist.m`
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/largedist.png)

## Working with RP (random projection) models
`randomprojection.m`
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/randomprojection.png)

## Specifying a custom list of correlations
`correlation.m`
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/correlations.png)

## Constructing composite models
`composite.m`
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/composite.png)

## Constructing and sampling a time-dependent model from high order Markov chains.
`timedep.m`
![](https://github.com/HussainAther/neuroscience/blob/master/maxent/schneidman/doc/timedep.png)
