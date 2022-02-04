# ContingencyTableAnalysis

Implementation of various models for contingency table analysis

---

## Overview

In contingency table analysis, we are generally interested in the independence of the row and column categories.
However, in a square contingency table, where the categories of the row and column variables are ordered and consist of the same classifications, the observed frequencies tend to be concentrated around the diagonal cells, and the interrelatedness between the two classifications is so strong that statistical independence is often not valid.
Therefore, instead of independence between categories, an association model that shows the structure of the odds ratio or a model with symmetry (or asymmetry) that shows the structure of the ratio to the probability of symmetrically located cells have been proposed.
Among the models with an association, there are the uniform association (U) model and the quasi-uniform association (QU) model proposed by Goodman (1979).
On the other hand, models with symmetry (or asymmetry) include the symmetry (S) model (Bowker, 1948), the quasi-symmetry (QS) model (Caussinus, 1965), and the k-th linear asymmetry (LSk) model (Tahata and Tomizawa, 2011) and so on.
The symmetry plus quasi-uniform association (SQU) model (Yamamoto and Tomizawa, 2010) and the k-th linear asymmetry plus quasi-uniform association (LSQUk) model (Fujisawa and Tahata, 2022) have been proposed as models with both properties.

In this project, we can get the output of applying different models to the data.
By comparing the values there, we get different results and interpretations.

---

## Documentation

There are two functions in models.R for analysis: `model` and `detail`.

### `model(freq, sort)`

This function displays the results of the analysis with various models.
This takes as its argument the table data `freq` that you want to analyze, and can optionally sort it in ascending order by adding the argument `sort` to any of the degrees of freedom, likelihood ratio chi-squared statistic, AIC, or p-value.
The sort argument can be given as string "df", "G^2", "AIC", or "Pr(>G^2)" as needed, respectively.

### `detail(model)`

The detail function returns detailed results for models of interest from the results of various models obtained by the model function.
This takes model as an argument.

---

## Operating environment

This was implemented using `R 4.0.5` with the `Rsolnp 1.16` and `stats 4.0.5` (especially, we use the `glm` function) packages.
