# cfisher: Fisher's exact test controlling for confounders

Fisher's exact test compares a binary exposure with a binary outcome and gives an exact
p-value, but it cannot adjust for confounders. In practical studies, confounders are
ubiquitous, and they can drive both the exposure and the outcome and create a false positive
association. The traditional logistic regression adjusts for confounders, but its p-values
are approximate and the test often fails when the exposure is rare or the data are separated.

`cfisher.test()` keeps the exact test and adds the control for confounders. It calculates
each subject's probability of being a case under the null hypothesis of no association,
using a logistic model of the outcome on the confounders without the exposure. It then
computes the exact null distribution of the number of exposed cases, and returns the
p-value, the conditional odds ratio estimate, and its confidence interval in the same
`htest` format as `fisher.test()`. So users who are familiar with `fisher.test()` can easily
adopt the function by including confounder data. Without confounders the results are those
of `fisher.test()`. With confounders the fitted model is plugged in as if known, so the
p-value and the interval are exact at the true effects of the confounders and approximate
at the fitted ones. The null model can be fitted once and passed to many tests, which suits
a scan of many exposures (for example, genetic variants) against one outcome.

## Installation

The package is on GitHub (repository `ZWuLab/cfisher`):

```r
# install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github("ZWuLab/cfisher", build_vignettes = TRUE)
```

The package contains compiled C++ code, so a C++ compiler is needed for this installation.
`build_vignettes = TRUE` installs the vignette as well, which needs `knitr` and `rmarkdown`.

## Example

A confounder `z` raises both the chance of exposure and the chance of being a case. The
exposure has no effect on the outcome.

```r
library(cfisher)
set.seed(2026)
n <- 600
z <- rnorm(n)                                # a confounder
x <- rbinom(n, 1, plogis(-1 + 1.2 * z))      # the confounder raises the chance of exposure
y <- rbinom(n, 1, plogis(-0.5 + z))          # ... and of being a case; x itself has no effect
fisher.test(table(x, y))                     # ignores z
#>
#> 	Fisher's Exact Test for Count Data
#>
#> data:  table(x, y)
#> p-value = 3.767e-05
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  1.444957 2.956233
#> sample estimates:
#> odds ratio
#>   2.064752
cfisher.test(x, y, z)                        # controls for z
#>
#> 	Fisher's exact test controlling for confounders; odds ratio: slopes
#> 	refitted under each hypothesized odds ratio
#>
#> data:  x and y, controlling for z
#> T = 108, exposed = 206, cases = 245, n = 600, p-value = 0.453
#> alternative hypothesis: true odds ratio is not equal to 1
#> 95 percent confidence interval:
#>  0.4909373 1.3388122
#> sample estimates:
#> conditional odds ratio
#>              0.8215806
```

At the 1% level, Fisher's exact test, which ignores `z`, reports an association (p-value
3.8e-05), which is wrong. The test that controls for `z` does not (p-value 0.45). The two
results print in the same layout, so their estimates and intervals can be compared as well.

## More hints

* `cfisher.test(x, y, z, conditional = FALSE)` uses the unconditional null distribution
  (Poisson-binomial) instead of the conditional one.
* `cfisher.test(x, y, fit = fit)` takes a fitted binomial `glm`, or a vector of null
  probabilities, in place of `z`. With a vector of probabilities the slopes cannot be
  refitted. The odds ratio then uses the exponential tilt of the distribution at those
  probabilities, which is exact only when they are the baseline probabilities under the
  alternative (otherwise not recommended).
* `conf.int = FALSE` omits both the odds ratio estimate and the interval, for speed in
  scans (`fisher.test()` still returns its estimate).
* `cfisher.null()` returns the null distribution of the number of exposed cases.
* The vignette `vignette("cfisher")` explains the test, how to read its output and how to
  scan many exposures.

## Reference

Wu, Z. and Contador, G. (2026). Fisher's exact test controlling for confounders. Draft
manuscript.

> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
