# Re-creating weak pairwise correlations

This repository contains the code and output in recreating the work of Schneidman et al., (2006). "Weak pairwise correlations imply strongly correlated network states in a neural population." Nature.

The scripts use code by the [maxent_toolbox](https://github.com/orimaoz/maxent_toolbox).

* Data in `data`.
* Code in `src`.
* Output in `doc`.

## Working with small distributions of neurons (exhaustively)
Used in `smalldist.m`.
![](https://raw.githubusercontent.com/HussainAther/neuroscience/master/maxent/schneidman/doc/smalldist.png)

## Working with larger distributions of neurons (MCMC)
`largedist.m`

## Working with RP (random projection) models
`randomprojection.m`

## Specifying a custom list of correlations
`correlation.m`

## Constructing composite models
`composite.m`

## Constructing and sampling a time-dependent model from high order Markov chains.
`timedep.m`
