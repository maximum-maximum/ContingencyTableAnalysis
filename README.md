# ContingencyTableAnalysis

Implementation of various models for contingency table analysis

## Overview

In contingency table analysis, we are generally interested in the independence of the row and column categories.
However, in a square contingency table, where the categories of the row and column variables are ordered and consist of the same classifications, the observed frequencies tend to be concentrated around the diagonal cells, and the interrelatedness between the two classifications is so strong that statistical independence is often not valid.
Therefore, instead of independence between categories, an association model that shows the structure of the odds ratio or a model with symmetry (or asymmetry) that shows the structure of the ratio to the probability of symmetrically located cells have been proposed.
Among the models with an association, there are the uniform association (U) model and the quasi-uniform association (QU) model proposed by Goodman (1979).
On the other hand, models with symmetry (or asymmetry) include the symmetry (S) model (Bowker, 1948), the quasi-symmetry (QS) model (Caussinus, 1965), and the k-th linear asymmetry (LSk) model (Tahata and Tomizawa, 2011) and so on.
The symmetry plus quasi-uniform association (SQU) model (Yamamoto and Tomizawa, 2010) and the k-th linear asymmetry plus quasi-uniform association (LSQUk) model (Fujisawa and Tahata, 2022) have been proposed as models with both properties.

In this project, we can get the output of applying different models to the data.
By comparing the values there, we get different results and interpretations.

## Documentation

There are two functions in models.R for analysis: `model` and `detail`.

### `model(freq, sort)`

This function displays the results of the analysis with various models.
This takes as its argument the table data `freq` that you want to analyze, and can optionally sort it in ascending order by adding the argument `sort` to any of the degrees of freedom, likelihood ratio chi-squared statistic, AIC, or p-value.
The sort argument can be given as string "df", "G^2", "AIC", or "Pr(>G^2)" as needed, respectively.

### `detail(model)`

The detail function returns detailed results for models of interest from the results of various models obtained by the model function.
This takes model as an argument.

## Example

For an example analysis, consider the following square contingency table.

| ＼  | Y_1 | Y_2 | Y_3 | Y_4 |
| --- | --- | --- | --- | --- |
| X_1 | 6   | 2   | 3   | 1   |
| X_2 | 9   | 4   | 2   | 1   |
| X_3 | 9   | 2   | 3   | 1   |
| X_4 | 12  | 1   | 2   | 1   |

First, we represent the data in the form of a vector and apply the `model` function.
At this time, we also sort the output results by AIC.

```
> freq <- c(6, 2, 3, 1, 9, 4, 2, 1, 9, 2, 3, 1, 12, 1, 2, 1)
> model(freq, "AIC")
```

Then, we get the following results.

```
   model      df   G^2   AIC Pr(>G^2)
6   LSI1      11  7.37 63.18   0.7680
8   LSI3  (I)  9  3.72 63.53   0.9287
7   LSI2      10  6.30 64.11   0.7893
9   LSU1      10  6.32 64.13   0.7877
11  LSU3  (U)  8  2.69 64.49   0.9524
10  LSU2       9  5.27 65.08   0.8101
12 LSQI1       7  3.62 67.42   0.8228
14 LSQI3 (QI)  5  0.77 68.57   0.9792
13 LSQI2       6  2.98 68.79   0.8114
24 Cov=0       1  0.98 68.87   0.3211
15 LSQU1       6  3.61 69.41   0.7299
17 LSQU3 (QU)  4  0.69 70.49   0.9529
18   LS1       5  2.97 70.78   0.7042
16 LSQU2       5  2.98 70.79   0.7031
19   LS2       4  2.33 72.13   0.6755
20   LS3 (QS)  3  0.46 72.26   0.9283
2     SU      11 22.53 78.34   0.0206   *
1     SI      12 26.98 80.79   0.0078  **
3    SQI       8 19.98 81.79   0.0104   *
22   ME2       2 17.08 82.97   0.0002 ***
21   ME3 (MH)  3 19.12 83.01   0.0003 ***
4    SQU       7 19.86 83.66   0.0059  **
23   ME1       1 16.43 84.32   0.0001 ***
5      S       6 19.27 85.07   0.0037  **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

From these results, we want to know more about the analysis of the LSI1 model, so we apply the `detail` function.

```
> detail("LSI1")
```

```
Model: LSI1

Call:
glm(formula = formula, family = poisson, data = list(freq))

Deviance Residuals:
     Min        1Q    Median        3Q       Max
-1.62470  -0.32497   0.03754   0.46459   0.83855

Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  0.03044    0.50187   0.061 0.951631
arrayLSI11   1.03211    0.29503   3.498 0.000468 ***
arrayLSI12   0.47071    0.32822   1.434 0.151540
arrayLSI13   0.41350    0.31340   1.319 0.187034
arrayLSI14  -0.39004    0.09517  -4.098 4.16e-05 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 41.6191  on 15  degrees of freedom
Residual deviance:  7.3743  on 11  degrees of freedom
AIC: 63.18

Number of Fisher Scoring iterations: 4

Data:
     [,1]        [,2]       [,3]       [,4]
[1,] "6 (8.12)"  "2 (3.14)" "3 (2.01)" "1 (0.9)"
[2,] "9 (6.84)"  "4 (2.64)" "2 (1.69)" "1 (0.76)"
[3,] "9 (9.55)"  "2 (3.69)" "3 (2.36)" "1 (1.06)"
[4,] "12 (9.32)" "1 (3.6)"  "2 (2.3)"  "1 (1.03)"
(The parenthesized values are the MLEs of expected frequencies under the selected model)
```

By doing so, we can get the detailed analysis results for each parameter and the MLE values.

## Operating environment

This was implemented using `R 4.0.5` with the `Rsolnp 1.16` and `stats 4.0.5` (especially, we use the `glm` function) packages.
