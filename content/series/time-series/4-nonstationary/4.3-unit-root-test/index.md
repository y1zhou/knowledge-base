---
title: "Unit Root Test"
date: 2020-10-07T16:42:23-04:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "time-series-unit-root-test"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 43 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

We've been applying different orders of differencing on simulated data by looking at the time series plot and ACF plot. In this chapter, we introduce the `Dickey-Fuller unit-root test`, or the DF test, that helps us check whether differencing is necessary or not. The test shouldn't be used as the single source of truth, but rather an addition to the existing plotting methods to decide whether to apply differencing.

## Dickey-Fuller test

If $X_t$ is either an AR(1) model or an AR(1) with some mean trend:

$$
X_t = \phi X_{t-1} + Z_t, \quad X_t = \delta + \phi X_{t-1} + Z_t, \quad X_t = \delta_0 + \delta_1 t + \phi X_{t-1} + Z_t,
$$

then the null hypothesis of the unit-root test is $H_0: \phi = 1$ (nonstationary)[^unit-root] and the alternative hypothesis is $H_1: \phi < 1$ (stationary). If we fail to reject $H_0$, then the series might be a random walk and the stochastic trend could be removed by differencing. If we reject the null, then the time series is stationary and no differencing is necessary.

[^unit-root]: We're testing whether $\phi=1$, hence the name unit root.

The test statistic is

$$
DF = t\text{-ratio} = \frac{\hat\phi - 1}{se(\hat\phi)}
$$

where $\hat\phi$ is the LSE. The p-value could be calculated by a nonstandard asymptotic distribution, which is repeatedly generating data (resampling / bootstrapping) from a random walk to get the distribution of $\hat\phi$. In R,

```r
set.seed(1)
n.sim <- 1000
sample.size <- 100
phi.dist <- numeric(n.sim)

for (i in 1:n.sim) {
  e <- rnorm(sample.size)
  y <- cumsum(e)  # White noise
  ylag <- c(NA, y[-sample.size])
  dy.fit <- lm(y ~ ylag - 1)
  phi.hat <- coef(dy.fit)
  phi.hat.se <- coef(summary(dy.fit))[2]
  phi.dist[i] <- (phi.hat - 1) / phi.hat.se
}

quantile(phi.dist, c(0.01, 0.025, 0.05, 0.1,
                     0.9, 0.95, 0.975, 0.99))
#        1%      2.5%        5%       10%       90%       95%     97.5%       99%
# -2.277308 -2.085253 -1.900637 -1.599258  0.863542  1.179292  1.426069  1.782729
```

## Augmented DF unit-root test

The DF test has a strong model assumption, and it's very rare for a real time series to be explained by a random walk. If we observe an increasing pattern in the time series plot, and a slowly decaying pattern in the ACF plot, but the DF test says no differencing is needed, then maybe the data is not generated from the specific model.

The `augmented Dickey-Fuller` (ADF) test relaxes the assumption to an AR(p) model:

$$
\nabla X_t = c_t + \beta X_{t-1} + \sum_{i=1}^{p-1} \phi_i \nabla X_{t-i} + Z_t
$$

The null hypothesis is $H_0: \beta_0$ (nonstationary), and the alternative hypothesis is $H_1: \beta < 0$. Similarly, if we fail to reject the null, then differencing might be necessary. If we reject $H_0$, then the series might be stationary or is an explosive AR[^adf-test-types].

[^adf-test-types]: There are two types of ADF tests, depending on the alternative hypotheses.

In the following section we're going to see how the ADF test works on different datasets.

```r
library(tseries)

set.seed(1)
x <- arima.sim(list(order = c(1, 0, 0), ar = 0.7), n = 100)
adf.test(x)

# 	Augmented Dickey-Fuller Test

# data:  x
# Dickey-Fuller = -4.375, Lag order = 4, p-value = 0.01
# alternative hypothesis: stationary

# Warning message:
# In adf.test(x) : p-value smaller than printed p-value
```

By default, the `alternative` parameter is set as "stationary" in the `adf.test()` function. With a p-value smaller than 0.01, we reject the null and conclude that the process is stationary. We could also check if the data is an explosive AR:

```r
adf.test(x, alternative = "explosive")

# 	Augmented Dickey-Fuller Test

# data:  x
# Dickey-Fuller = -4.375, Lag order = 4, p-value = 0.99
# alternative hypothesis: explosive

# Warning message:
# In adf.test(x, alternative = "explosive") :
#   p-value smaller than printed p-value
```

And clearly it's not. What if we use this on an actual explosive AR series?

```r
data("explode.s", package = "TSA")
adf.test(explode.s, alternative = "explosive")

# 	Augmented Dickey-Fuller Test

# data:  explode.s
# Dickey-Fuller = 6.1576, Lag order = 1, p-value = 0.01
# alternative hypothesis: explosive

# Warning message:
# In adf.test(explode.s, alternative = "explosive") :
  # p-value greater than printed p-value
```

With a p-value smaller than 0.01, we reject the null and conclude that we have a nonstationary explosive AR series.

## Over-differencing

Unnecessary differencing can cause artificial autocorrelation and make the model fitting process complicated. For example, suppose $X_t$ is a random walk, then

$$
\nabla X_t = X_t - X_{t-1} = Z_t
$$

But what happens after second order differencing?

$$
\nabla^2 X_t = \nabla(Z_t) = Z_t - Z_{t-1}
$$

With first order differencing we can simply fit a white noise, but after another round differencing an MA(1) model would be needed. Always keep the `principle of parsimony` in mind!

Another example with $X_t \sim WN$. With some manipulation,

$$
\begin{gathered}
  X_t = Z_t, X_{t-1} = Z_{t-1} \\\\
  -0.5 X_{t-1} + 0.5 Z_{t-1} = 0 \\\\
  \Rightarrow X_t = -0.5 X_{t-1} + Z_t + 0.5 Z_{t-1}
\end{gathered}
$$

The WN could be modeled as an ARMA(1, 1). If we have a choice between the two, WN is preferred by the principle of parsimony.
