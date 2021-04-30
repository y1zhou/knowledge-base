---
title: "A Bayesian Perspective on Missing Data Imputation"
date: 2021-04-26T10:20:00-05:00
summary: 'This lecture discusses some approaches to handling missing data, primarily when missingness occurs completely randomly. We discuss a procedure, MICE, which uses Gibbs sampling to create multiple "copies" of filled-in datasets.' # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-missing-data-imputation"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 130 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Missing piece in a jigsaw." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "B-x4VaIriRc" # Unsplash ID of the picture
---



Missing data is a quite prevalent issue among many applied data sets. This lecture is a light introduction to one way we can perform missing data imputation from a Bayesian perspective. For more information, the book *Statistical Analysis with Missing Data* by Roderick Little and Donald Rubin is a great resource.


## Diabetes example

The motivating example comes from Hoff's book, page 115. Consider two measurements on 200 Arizona women with Pima Indian heritage:

- `skin`, the skin fold thickness, and
- `bmi`, the body mass index.

Our goal is to study the relationship between these two variables, accounting for the missing data present in the table.


```r
library(tidyverse)
library(ggpubr)
library(rjags)
```


```r
diab <- data.table::fread(
  "http://www2.stat.duke.edu/~pdh10/FCBS/Exercises/diabetes_200_miss.dat",
  header = T
) %>%
  select(skin, bmi)

diab %>%
  head() %>%
  knitr::kable()
```



| skin|  bmi|
|----:|----:|
|   28| 30.2|
|   33|   NA|
|   NA| 35.8|
|   43| 47.9|
|   NA|   NA|
|   27|   NA|


### Exploratory data analysis

As shown in the table above, there are obviously missing values coded as NA's in our data set[^missing-data-code].

[^missing-data-code]: This is **not** the only way missing data is coded. For example, in some software -999 could be used to represent a missing value since it's "obviously" unrealistic.

We can use the `summary()` function in R to get the 5-number summary, mean, and number of missing values for the columns.


|        | skin|   bmi|
|:-------|----:|-----:|
|Min.    |  7.0| 18.20|
|1st Qu. | 20.5| 27.60|
|Median  | 29.0| 32.80|
|Mean    | 29.1| 32.21|
|3rd Qu. | 36.0| 36.48|
|Max.    | 99.0| 47.90|
|NA's    | 25.0| 22.00|

Histograms[^hist-code] are helpful in determining the distributions of the variables. Both are roughly symmetric, and there's one extreme observation for skin fold thickness.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/diab-hist-1.svg" alt="Observed BMI measurements and skin fold thicknesses." width="576" />
<p class="caption">Figure 1: Observed BMI measurements and skin fold thicknesses.</p>
</div>

[^hist-code]: R code for generating Figure <a href="#fig:diab-hist">1</a>:

    
    ```r
    diab %>%
      pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
      mutate(
        Label = case_when(
          Variable == "skin" ~ "Observed skin fold thickness",
          T ~ "Observed BMI measurements"
        )
      ) %>%
      gghistogram(x = "Value", bins = 10, xlab = "") %>%
      facet(facet.by = "Label", scales = "free")
    ```

We can also use a scatter plot to show the relationship between the two variables.

<img src="{{< blogdown/postref >}}index_files/figure-html/diab-scatter-1.svg" width="480" />

We observe a relatively strong positive relationship except for the one outlier, and we might want to quantify this in terms of how correlated the two variables are.

Finally, when we have missing values in our data set, it can be helpful to visualize the **pattern of missingness**. The `aggr()` function from the `VIM` R package shows proportions of missing values for each variable, as well as how frequently missingness occurrs for various variable combinations.


```r
VIM::aggr(diab)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/diab-missing-pattern-1.svg" width="576" />

There's about 13% missingness for the `skin` variable, and about 11% for the `bmi` variable. The blue color on the right indicates presence of data, whereas the red indicates missing values. We have both values for the large majority of the data, but in a few entries we are missing one or both of the values.

Different strategies might be applied for different types of missingness. For example, if we *always* have values in one variable but have missing entries in another variable, or if the two variables are always missing together.


## Bayesian models

Our two variables appear correlated with each other, so we specify a `bivariate normal` data distribution for the **sampling model**:

<div>
$$
  \boldsymbol{Y}_i = \begin{bmatrix} Y_{i1} \\ Y_{i2} \end{bmatrix} \mid \theta, \Sigma
  \overset{iid}{\sim} 
  \text{MVN} \left(
    \textcolor{crimson}{\begin{bmatrix} \theta_1 \\ \theta_2 \end{bmatrix}},
    \textcolor{orange}{\begin{bmatrix} \Sigma_{11} & \Sigma_{12} \\ \Sigma_{21} & \Sigma_{22} \end{bmatrix}}
  \right)
$$
</div>

where `\(Y_{i1}\)` is the `\(i^{th}\)` person's skin, and `\(Y_{i2}\)` is the `\(i^{th}\)` person's BMI. The two parameters of the MVN distribution are its <span style="color: crimson">mean vector</span> and <span style="color: orange">covariance matrix</span>. The multivariate normal is just an extension of the normal distribution that accounts for correlation between components of the random variables[^separate-normals].

[^separate-normals]: If we don't assume relationship between the two variables, then we could easily specify separate normal distributions for each.

For the **prior model**, we have essentially six parameters. `\(\theta_1\)` and `\(\theta_2\)` are means, so we specify a normal distribution with a large variance:

<div>
$$
  \theta_1, \theta_2 \overset{iid}{\sim} \mathcal{N}(0, 10^{10})
$$
</div>

For the covariance matrix of the multivariate normal distribution, there's a conjugate prior called the `inverse Wishart distribution`, which is a distribution over matrices that has two parameters.

**Posterior distributions** for our parameters are given by Bayes' Theorem:

`\begin{aligned}
  p(\theta, \Sigma \mid y) &\propto p(y \mid \theta, \Sigma) p(\theta, \Sigma) \\
  &\propto \left[ \prod_{i=1}^n p(y_i \mid \theta, \Sigma) \right] p(\theta) p(\Sigma)
\end{aligned}`

A problem arises when we evaluate the sampling model likelihood, since this depends on **all** of the data, but we have missing values. For example, if we were missing `\(Y_{31}\)`, then we can't evaluate `\(p(y_3 \mid \theta, \Sigma)\)`.


## Simple approaches

There are many ways to handle missing data. The simplest approach is to perform a `complete case analysis`, where we delete entries with a missing value for at least one variable. If the data are missing completely at random (`MCAR`), then we get unbiased parameter estimate. However:

- Credible intervals will usually be too wide, since the sample size has effectively been reduced (drastically).
- If missing values occur for other reasons, this yields biased estimates.

The next step up is instead of dropping all the missing entries, we replace missing values of a feature with the mean or median of available values for the feature. This is called `single mean/median imputation`, and is often recommended in many machine learning tasks if the data are MCAR as it can still lead to unbiased parameter estimates.

There are also some drawbacks with this method, especially if you're interested in studying the variability. Single mean/median imputation artificially reduces the variance of features, which results in credible intervals which are too narrow. As shown in Figure <a href="#fig:bmi-imputed">2</a>, the peak at the mean is much higher relative to the surrounding values. The mean of the distribution can still be estimated appropriately, but the amount of spread is reduced.


```r
diab %>%
  select(bmi) %>%
  mutate(Imputed = bmi) %>%
  replace_na(list(Imputed = mean(.$Imputed, na.rm = T))) %>%
  select(`No imputation` = bmi, `With mean imputation` = Imputed) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "BMI") %>%
  gghistogram(x = "BMI", bins = 10, facet.by = "Variable")
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/bmi-imputed-1.svg" alt="Mean imputation with the BMI feature." width="576" />
<p class="caption">Figure 2: Mean imputation with the BMI feature.</p>
</div>



Another drawback is it doesn't account for relationship *between variables*, thus reduces correlation. For example, the correlation between the two variables is 0.665 if we only look at the complete cases. If we use single mean imputation, the correlation decreases to 0.589.

Note that if you don't care about assessing variability in your estimates, which is often the case for classification tasks, then mean imputation can work quite well.

## A Bayesian approach

A popular approach to address the issues of single mean imputation is called multiple imputation by chained equations (`MICE`), and proceeds as follows.

1. Specify separate **imputation models** for each feature with missing values, which are conditional on all other features. This allows us to exploit relationship between variables.
2. Use Gibbs sampling with the imputation models to generate multiple imputed data sets. This has the effect of allowing us to more appropriately assess uncertainty in the missing values.
3. Analyze each imputed data set separately, i.e. by sampling from the posterior using JAGS.
4. Combine the results across imputed data sets at the end.

Recall that in the [Gibbs sampler]({{< ref "/series/bayesian-stat/7-two-param-and-gibbs-sampling/index.md#gibbs-sampler" >}}) we generate values at step `\(s+1\)` using the full conditional distributions given values of other parameters at step `\(s\)`. In our case since missing values are unknown (just like parameters), we can replace the parameters in the Gibbs sampler with our missing values.

Instead of deriving full conditional distributions, MICE takes a shortcut for missing data by directly specifying full conditional distributions. One hopes that there is a proper posterior distribution which would yield these full conditional distributions.

A common choice is to specify a `Bayesian linear regression` model, in our to exploit relationships between variables. Let `\(Y_{i1}\)` denote the skin variable and `\(Y_{i2}\)` the BMI:

`\begin{gathered}
  Y_{i1} \mid Y_{i2}, \beta^{(2)}, {\sigma^2}^{(2)} \sim \mathcal{N}(\beta_1^{(2)} + \beta_2^{(2)} Y_{i2}, {\sigma^2}^{(2)}) \\
  Y_{i2} \mid Y_{i1}, \beta^{(1)}, {\sigma^2}^{(1)} \sim \mathcal{N}(\beta_1^{(1)} + \beta_2^{(1)} Y_{i1}, {\sigma^2}^{(1)})
\end{gathered}`

where the two mean terms are linear regression of BMI on skin and vice versa. Then we specify appropriate priors on `\(\beta^{(1)}\)`, `\(\beta^{(2)}\)`, `\({\sigma^2}^{(1)}\)`, and `\({\sigma^2}^{(2)}\)`. This gives us a Bayesian model to estimate the missing skin values `\(Y_{i1}\)` and missing BMI values `\(Y_{i2}\)`.

What's a bit different from before is that we are only generating a small number of samples, because we are filling in missing data instead of parameter values. This is usually only run for 5-10 iterations.

### Implementation in R

This imputation procedure can be done automatically in R using the `mice()` function in the `mice` package. Below we run `\(D=5\)` iterations to obtain `\(D\)` separate, imputed data sets:


```r
set.seed(42)
D <- 5
diab_mice <- mice::mice(
  diab,            # dataset with missing values
  m = D,           # number of imputed datasets
  method = "norm"  # use Bayesian linear regression discussed above
)

Y_imp <- map(seq(D), ~mice::complete(diab_mice, .x))
```

The `Y_imp` variable is a list of five data frames, each with missing values filled in by the corresponding values in `diab_mice$imp`.

Just like with other MCMC methods, we can assess convergence by looking at trace plots. With only five iterations, we want to see them randomly dispersed among each other, e.g. it would be undesirable if the red curve was always above the other curves.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/mice-trace-plots-1.svg" alt="Mean and standard deviation of imputed values across the five generated data sets." width="576" />
<p class="caption">Figure 3: Mean and standard deviation of imputed values across the five generated data sets.</p>
</div>

### Drawing posterior samples

We can then treat each of the five imputed data sets **separately**, and use JAGS[^diab-jags] to generate separate posterior samples for the parameters of interest. The table below has 95% credible intervals for the specified parameters:



[^diab-jags]: R code for running JAGS:

    
    ```r
    n <- nrow(Y_imp[[1]])  # sample size
    
    # Specify parameters, initial values, and JAGS settings
    parameters <- c("theta", "sig21", "sig22", "rho")
    
    initValues <- list(
      "theta" = unname(colMeans(diab, na.rm = T)),
      "tau" = matrix(c(1, 0, 0, 1), nrow = 2)
    )
    
    # JAGS settings
    adaptSteps <- 10000     # number of steps to "tune" the samplers
    burnInSteps <- 10000    # number of steps to "burn-in" the samplers
    nChains <- 3            # number of chains to run
    numSavedSteps <- 10000  # total number of steps in chains to save
    thinSteps <- 100        # number of steps to "thin" (1 = keep every step)
    nIter <- ceiling((numSavedSteps*thinSteps)/nChains) 	# steps per chain
    
    mcmcChain <- NULL
    for (d in seq(D)) {
      dataList <- list(
        "n" = n,
        "Y_imp" = as.matrix(Y_imp[[d]]),
        "I" = matrix(c(1, 0, 0, 1), nrow = 2)
      )
      
      m <- textConnection("
    model {
      # Sampling model - for multivariate normal, second argument is the
      # precision matrix, which is the inverse of the covariance matrix
      for (i in 1:n) {
      	Y_imp[i,1:2] ~ dmnorm(theta[], tau[, ])
      }
    
      # Prior model
      theta[1] ~ dnorm(0, 1e-10)
      theta[2] ~ dnorm(0, 1e-10)
      tau[1:2, 1:2] ~ dwish(I, 3)
    
      # Transformation
      # Obtain covariance matrix from precision matrix
      sigma[1:2, 1:2] <- inverse(tau[, ])
    
      # Extract variances and correlation
      sig21 <- sigma[1, 1]
      sig22 <- sigma[2, 2]
      rho <- sigma[1, 2] / (sqrt(sig21)*sqrt(sig22))
    }")
      jagsModel <- jags.model(m, 
                              data = dataList, 
                              inits = initValues, 
                              n.chains = nChains, 
                              n.adapt = adaptSteps)
      close(m)
      
      if (burnInSteps > 0) {
        update(jagsModel, n.iter = burnInSteps)
      }
      codaSamples <- coda.samples(jagsModel, 
                                  variable.names = parameters, 
                                  n.iter = nIter, 
                                  thin = thinSteps)
      mcmcChain[[d]] <- as.matrix(codaSamples)
    }
    ```


|$\theta_1$     |$\theta_2$     |$\rho$       |
|:--------------|:--------------|:------------|
|(27.59, 31.02) |(31.24, 32.93) |(0.57, 0.73) |
|(27.89, 31.16) |(31.17, 32.86) |(0.58, 0.74) |
|(27.48, 30.71) |(31.25, 32.93) |(0.58, 0.74) |
|(27.84, 31.23) |(31.23, 32.91) |(0.61, 0.75) |
|(27.67, 30.88) |(31.62, 33.34) |(0.55, 0.71) |

These values are fairly similar, mainly because the sample size is large, and there is a relatively small proportion of missing values.


### Combining results

Generally, we do not care about posterior inference of parameters for each imputed data set. Instead, we would rather **pool** our results together to account for imputation variability.

In frequentist inference, this is done using `Rubin's rules`; in Bayesian inference, we just aggregate all of our posterior samples together across the imputed data sets. The credible intervals after aggregating posterior samples are given below.


|           |95\% CI        |
|:----------|:--------------|
|$\theta_1$ |(27.67, 31.04) |
|$\theta_2$ |(31.26, 33.07) |
|$\rho$     |(0.57, 0.74)   |

From the 95% CI of `\(\rho\)`, we conclude there is fairly strong positive correlation between skin and BMI.


### Concluding remarks

For many inferential tasks, simple imputation methods can perform well! Advanced methods are necessary if one is focused on quantifying uncertainty in estimates correctly, particularly for data sets with low sample sizes or a large amount of missingness.

Missing data is very application-specific, so you should really think about what the variables represent, i.e. use domain knowledge and explore the data before considering how to impute. Missing data often does **not** occur completely randomly, e.g.:

- Non-response to survey questions -- it might be a sensitive question, or maybe the survey is too long.
- Patient dropout in a longitudinal study -- maybe the patients died or moved to somewhere else.

Finally, **never, ever replace all of the missing values with zeros (or other arbitrary values)**!
