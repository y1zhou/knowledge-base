---
title: "Bootstrap"
slug: "nonparametric-methods-bootstrap"
categories:
  - Nonparametric Methods
tags:
  - Nonparametric Methods
  - Statistics
  - Sampling

summary: "The procedure and applications of the nonparametric bootstrap."
date: 2019-05-06T10:46:18-04:00
toc: true
weight: 51
---

Many of the modern nonparametric statistics are highly computational - mostly been developed in the last 30 years or so.

Bootstrapping is a widely used computational method, mainly because often the theory is intractable - you can't easily work out the properties of the estimator, or the theory is difficult to apply / work with - checking the assumptions itself can be difficult. It allows us to work with

- (almost) arbitrary distributions at the population level
- more general characteristics of the population, i.e. not just the mean

## Introduction
There are two main versions of the bootstrap: the parametric bootstrap and the `non-parametric bootstrap`.

A nice thing about the nonparametric bootstrap is that it's very close to being **assumption free**. One key assumption is that the sample (empirical) CDF is a good approximation to the true population CDF for samples that are not too small. That is, $S(x)$, the ECDF, should reflect the main characteristics of $F(x)$, the theoretical population CDF. We don't care about $F(x)$ - we want to use $S(x)$ as a **proxy** for $F(x)$.

The idea is based on `resampling` of data. In the real world, we have the population $F(x)$ and some parameter $\theta$. Repeated sampling is usually not feasible. In the bootstrap world, we have a sample $S(x)$ and an estimated parameter $\hat{\theta}$. Repeated sampling is totally feasible with enough computing power to do it.

Note that permutation testing is also a type of resampling. So what's the main difference? With `permutation`, we're resampling without replacement and leads to exact tests. With `bootstrap`, we resample *with replacement* from our original sample, and that leads to approximate results. Bootstrapping leads to a more general approach.

## Procedure
Suppose we have a sample of size $n$: $x_1, x_2, \cdots, x_n$ from some unknown population. A `bootstrap sample` would be a random sample of size $n$ sampled from $x_1, x_2, \cdots, x_n$ with replacement. Some $x_i$ will appear multiple times in a given bootstrap sample, and others not at all.

In principle, we could enumerate all possible bootstrap resamples; in practice, this number grows very quickly with $n$, e.g. when $n=10$ there are more than 90,000 possible different bootstrap samples.

In real applications, we'll always take a sample of the possible bootstrap resamples - call this $B$. Often, $B \approx 500$ will suffice, or $B \approx 1000$ but we rarely need more than this.

**Standard notation**: a typical bootstrap resample will be denoted $x_1^\ast, x_2^\ast, \cdots, x_n^\ast$, where each $x_i^\ast$ is equal to one of the original $x_i$. If $B$ bootstrap resamples are generated, the $b^{th}$ is denoted by $x_1^{\ast b}, x_2^{\ast b}, \cdots, x_n^{\ast b}$, $b = 1, \cdots, B$.

For each bootstrap sample, we compute the statistic of interest (getting $B$ values, one for each sample). We build up from these an *approximate, empirical* distribution. We use this to learn about the population quantity of interest, e.g. estimate the population parameter by the **average** of the $B$ bootstrap quantities.

Suppose $\theta$ is the population parameter of interest. Let $S(x^\ast)$ be an estimate of $\theta$ from the bootstrap sample $x^\ast$.  The **average** of the $S(x^\ast)$ values is an estimate for $\theta$.
$$
\hat{\theta} = S(\cdot^*) = \frac{1}{b}\sum_{b=1}^B{S(x^{\ast b})}
$$

>  The book call it $S(\cdot^\ast)$, which is non-typical; it's usually denoted $\hat{\theta}$ or $\hat{\theta^\ast}$.

As $B$ increases, $S(\cdot^\ast)$ tends to the mean computed for the **true** `bootstrap sampling distribution`, which is based on all possible configurations. Even for $B \approx 500$ or $1000$, the estimator will be pretty good in general.

Similarly, the estimator of the true bootstrap standard error is

$$
se_B(S) = \sqrt{\frac{1}{B-1} \sum_{b=1}^B{[S(x^{\ast b}) - S(\cdot^\ast)]^2}}
$$


where $S(x^{\ast b})$ is the sample s.d. of bootstrap sample values. The easiest way to do this in R is to use the `sample()` function with `replace = T`, e.g. `sample(x, size = n, replace = T)` where `x` is the data with `n` observations.

## Example in R
Suppose we have a small dataset with $n=20$ observations from some unknown distribution. We consider bootstrap for the mean and the median.

$$
\begin{gather*}
  0.049584203 && 0.165223813 && 0.070872759 && 0.009547975 \\\\
  0.118294157 && 0.365906785 && 0.078496145 && 0.102532392 \\\\
  0.297899453 && 0.178715619 && 0.336178899 && 0.602171107 \\\\
  0.024036640 && 0.014285340 && 0.313490262 && 0.077453934 \\\\
  0.118797809 && 0.155527571 && 0.311395078 && 0.092584865
\end{gather*}
$$


In the sample, we see $\bar{x} = 0.174$ and a median of $0.119$, which indicates skewness. A boxplot would also show this. We should be suspicious of a normal assumption, but the sample is small so it's hard to know. Maybe we don't want to use the mean (and its standard error) as our basis for inference, especially in light of the skewness. The median could be a better alternative, but we know little about the sampling distribution of the median (standard error, etc.). Some of these problems can be addressed via the bootstrap.

```r
B <- 500
boot.mean <- vector("numeric", length = B)
boot.median <- vector("numeric", length = B)
for (i in seq(B)) {
    boot.mean[i] <- mean(sample(x, size = length(x), replace = T))
    boot.median[i] <- median(sample(x, size = length(x), replace = T))
}
```

The mean and median of the `boot.mean` (0.1747501 and 0.1731665, respectively) are both pretty similar to the data mean, 0.1741497. The bootstrap standard deviation `sd(boot.mean)` 0.03444655 matches pretty well with the standard error for the mean `sd(x) / sqrt(length(x))` 0.03395322. The bootstrap mean behaves very much like what we would "traditionally" do for a data set like this, i.e. compute the sample mean and the standard error.

What we couldn't do so easily is to look at the median of the sample. We can use the sample median as an estimator for the population median, but we don't have an easy estimator for the variability of the estimator! With bootstrap, it's quite easy to explore.

The `mean(boot.median)` 0.1290285 and `median(boot.median)` 0.118546 are both pretty close to the sample median 0.118546. However, we can also easily get, for example, `sd(boot.median)` to be 0.03989125, or the distribution of the minimum, etc.

The data were generated from an exponential with rate $\lambda = 4$[^seed]. The true mean($\theta$) is $\lambda^{-1} = 0.25$, and the true median is $\ln{2} / \lambda = 0.174$[^exp-median], so both are within the bootstrap mean and median results. We've underestimated both slightly, because the sample values of the mean and median were on the small side, and that introduces **bias** in bootstrap estimates.

## Regression analysis with bootstrap
Suppose we have multivariate data $(x_1, y_1)$, $(x_2, y_2)$, $\cdots$, $(x_n, y_n)$ where $x_i$ is the explanatory variable and $y_i$ is the response variable. In general, there are different types of bootstrap schemes that can be developed. These schemes depend on the types of assumptions you're willing to make about how the data might have been generated.

We'll consider simple linear regression here, but the same general ideas hold for multiple regression and generalized linear models, etc. Our model is
$$
Y_i = \beta_0 + \beta_1X_i + \epsilon_i \quad i = 1\cdots n
$$


where $\epsilon_i$ are uncorrelated, mean $0$ and constant variance $\sigma^2$. In general, the $x_i$ may be (i) fixed by design, (ii) randomly sampled, or (iii) simply observed. We're going to proceed as if they were fixed by design (it doesn't make that much difference).

The two common approaches for the bootstrap are `model-based resampling` (or `resampling residuals`) and `case-based resampling`. Note that there are different versions of these as well. Both methods are implemented in the R package `car::Boot()`.

### Model-based resampling
1. Set the $x^\ast$ as the original $x$ values.
2. Fit the `OLS` (ordinary least squares) to the linear model and find `residuals` (some version of it).
3. Resample the residuals, call them $\epsilon^\ast$.
4. Define

$$
y_i^\ast = \hat{\beta_0} + \hat{\beta_1}x_i + \epsilon^\ast
$$

where $\hat{\beta_0}$ and $\hat{\beta_1}$ are OLS estimators from original $(x_i, y_i)$ pairs; $x_i$ are the original data; $\epsilon^\ast$ are the resampled residuals, and $y^\ast$ are the new responses.

5. Fit OLS regression to $(x_1, y_1^\ast), \cdots, (x_n, y_n^\ast)$ to get the intercept estimate $\hat{\beta_0^{\ast b}}$, slope estimate $\hat{\beta_1^{\ast b}}$, and variance estimate ${S^2}^{\ast b}$ where $b = 1 \cdots B$.

This scheme doesn't trust the model that much, and has the advantage that it retains the information in the explanatory variables. There's arguements about which residuals to use (raw or studentized), but it doesn't make much of a difference in practice.

### Case-based resampling
1. Resample the data pairs $(x_i, y_i)$ directly, keeping the pairs together, with replacement. In R (or any other language), this can be done by sampling the indices.
2. FIt the OLS regression model and get $\hat{\beta_0^{\ast b}}, \hat{\beta_1{\ast b}}, {S^2}^{\ast b}$, $b = 1 \cdots B$.

This method assumes nothing about the shape of the regression function or the distribution of the residuals. All it assumes is that the observations are independent, which is an assumption we make all the time anyways. 

### Confidence intervals
In either case, we get an empirical bootstrap distribution of $\hat{\beta_0^\ast}, \hat{\beta_1^\ast}, {S^2}^\ast$ values. We could use these to get confidence intervals for the parameters of interest. This is a widely-useful application of the bootstrap in general - CIs when the theoretical derivation is difficult.

We'll have $B$ "versions" of the statistic of interest: $\hat{\theta}_1^\ast, \hat{\theta}_2^\ast, \cdots, \hat{\theta}_B^\ast$. This results in an `empirical distribution`. There are two simple, rather intuitive approaches to get CIs for $\theta$. One caveat is that it may not work well in practice - you may need to use "corrections" or improvements.

Recall that for a sample of size $n$ from a normal distribution with variance $\sigma^2$, the $95\\%$ CI for the mean, $\mu$, is given by $\bar{x} \pm 1.96\frac{\sigma}{\sqrt{n}}$. If $\sigma^2$ is unknown, which is typical, we replace $\sigma$ with $S$, the sample standard deviation, and $1.96$ by the appropriate quantile for $t_{n-1}$ (say $2.5\\%$ on both tails). We also use this when data are not necessarily normal, according to the central limit theorem.

The first procedure is to **mimic the "normal-type" confidence intervals**. We can do something similar using bootstrap-based estimates instead:

$$
\hat{\theta}^\ast \pm 2 \cdot se^\ast(\hat{\theta}^\ast)
$$


where $\hat{\theta}^\ast$ is the bootstrap estimate of $\theta$ (mean of bootstrap values), and $se^\ast(\hat{\theta}^\ast)$ is the bootstrap standard error estimate (standard deviation of bootstrap values). This forces the CI to be **symmetric** about $\hat{\theta}^\ast$.

The second approach is to **use the empirical bootstrap distribution**. Suppose we have a histogram of $\hat{\theta}_1^\ast, \hat{\theta}_2^\ast, \cdots, \hat{\theta}_B^\ast$. We would take the $2.5\\%$ of bootstrap values in the lower and upper tails. This does **not** enforce symmetry around $\hat{\theta}^\ast$. Here's an example in R for the second approach:

```r
set.seed(42)
x <- c(1, 2, 3, 5, 8)
y <- c(3, 7, 9, 7, 12)

res <- vector("numeric", 1000)
for (i in seq(1000)) {
  xb <- sample(x, size = length(x), replace = T)
  yb <- sample(y, size = length(y), replace = T)
  res[i] <- mean(xb) / mean(yb)
}
quantile(res, c(0.025, 0.975))
#      2.5%     97.5% 
# 0.2444444 0.9476842 
```

## Possible bias in bootstrap
In the **parametric bootstrap**, we assume that the data come from some known distribution but with *unknown* parameters. The idea of the parametric bootstrap is to use the data to estimate the parameters, and then draw resamples from the estimated distribution. For example, we draw samples from $N(\bar{x}, S^2)$ or $Pois(\bar{x})$. This method is not used much in practice because of the strong distribution assumption.

We can see that going from the parametric bootstrap to model-based resampling to case-based resampling, we're making fewer and weaker assumptions. In other words, we're getting more and more "secure".

So why don't we *always* use the safest bootstrap? In resampling cases, there's the usual **bias-variance tradeoff**. Generally speaking, if we compare the confidence intervals on the same data for the same parameters using the three methods, the parametric bootstrap would give the narrowest intervals, and the case-based resampling would yield the loosest bounds. If we're very certain (which almost never happens) that the overall model is correct, then the parametric bootstrap would work best.

### Jackknife

There are various procedures for dealing with the possible bias of a statistic, and `Jackknife` is one of them. We can use this to estimate bias and correct for it.

Suppose the sample is $x_1, x_2, \cdots, x_n$ and we're interested in some population characteristic $\theta$. We estimate it based on the entire sample by $\hat{\theta}$. We get $n$ **additional** estimates of $\theta$ by replacing the original sample by a set of subsamples, leaving out one observation in turn each time. The $i^{th}$ jackknife sample of size $n-1$ is given by:

$$
x_1, x_2, \cdots, x_{i-1}, x_{i+1}, \cdots, x_n \quad i=1, \cdots, n
$$


On each subsample, we estimate $\theta$ by $\hat{\theta}\_{(i)}$ (taking out the $i^{th}$ observation) and get $\hat{\theta}\_{(1)}$, $\hat{\theta}\_{(2)}$, $\cdots$, $\hat{\theta}\_{(n)}$. The **mean** of jackknife estimates is
$$
\hat{\theta}\_{(\cdot)} = \frac{1}{n}\sum\_{i=1}^n{\hat{\theta}\_{(i)}}
$$
and the Jackknife estimator is given as

$$
\hat{\theta}^J = n\hat{\theta} - (n-1)\hat{\theta}_{(\cdot)}
$$


where $n\hat{\theta}$ is based on the whole sample, and $(n-1)\hat{\theta}_{(\cdot)}$ is based on $n$ Jackknife estimates.

The `Jackknife estimator of bias` is simply defined as:

$$
\text{bias}(\hat{\theta}) = (n-1)\left(\hat{\theta}_{(\cdot)} - \hat{\theta}\right)
$$
and the bias-corrected jackknife estimate of $\theta$ is given by
$$
\hat\theta\_{jack} = \hat\theta^J - \text{bias}(\hat\theta) = n\hat\theta - (n-1)\hat\theta\_{(\cdot)}
$$

### Remarks
1. There are also purely bootstrap-based ways of dealing with bias, e.g. `Double bootstrap` (bootstrapping the bootstrap samples), or jackknife after bootstrap, etc.
2. The bootstrap is based on the "plug-in" principle, and the jackknife is based on the "leave one out" principle.
3. The Jackknife was developed before bootstrapping, and the latter is often considered as a superior technique. A comparison of the two methods is given [here](https://en.wikipedia.org/wiki/Resampling_%28statistics%29#Comparison_of_bootstrap_and_jackknife).

[^seed]:
    Generate the data with the following code:
    ```r
    set.seed(42)
    x <- rexp(20, 4)
    ```

[^exp-median]: Get this by solving $\int_0^m f(x)dx = 0.5$ where $f(x)$ is the PDF for the exponential distribution.

