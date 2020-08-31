---
title: "ARMA Models"
date: 2020-08-29T22:45:28-04:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "time-series-arma"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 20 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

In this chapter, we introduce some of the key properties of an important class of time series known as ARMA (autoregressive moving average) series.

## Autoregressive series

Let's consider a class of time series known as `autoregressive` (AR) series. The concept of autoregression related past observations of the series to current/future observations of the series using a regression-type model. We'll simulate 100 values from the following AR(1) model using R[^ar-r-code]:

[^ar-r-code]:
    For definition of the `plot_time_series` function, see [this post]({{< relref "1-introduction/index.md#fn:2" >}}).

    ```r
    gen_ar <- function(n, phi, sigma) {
    x <- rnorm(n, mean = 0, sd = sigma)
    x0 <- 0
    res <- NULL
    for (i in seq(n)) {
        x0 <- phi * x0 + x[i]
        res <- c(res, x0)
    }
    res
    }

    set.seed(1)
    x <- gen_ar(200, 0.7, 1)
    plot_time_series(x, x_type = "AR(1)")
    ```

$$
X_t = 0.7 X_{t-1} + Z_t,
$$

where $Z_t$ values are generated from the standard normal distribution. This is an AR model of order 1 because the current value depends on the value before that, called the `lag 1 value`. The coefficient $0.7$ tells us how strongly the time series depends on the previous value.

{{< figure src="ar_series.png" caption="Simulated AR data and its ACF and PACF." numbered="true" >}}

In the time series plot, focusing on timepoint 100 and onwards[^time-series-stabilization], notice how the correlation between adjacent values is high, i.e. neighboring values are close and not jumping around. This is not surprising as we generated the data by setting $X_t$ to be approximately $0.7 X_{t-1}$.

[^time-series-stabilization]: With simulated data, it's common to throw the first couple data points away because the time series has to be stabilized.

Also, there seems to be a trend of going up and down, and infrequently crossing the theoretical mean of zero (unlike white noise). This trend is due to the strong autocorrelation of neighboring values of the series.

The ACF plot gives the correlations between $X_t$ and $X_{t-1}$, $X_t$ and $X_{t-2}$ etc. Note how the correlation values decrease as the lag increases. Unlike the random walk, this is an `exponential decay` - the ACF declines very rapidly after some lags and essentially taper off to zero. {{<hl>}} Exponentially declining autocorrelations is a characteristic of autoregressive series.{{</hl>}}

The PACF spikes at lag 1 and is essentially zero for other lags. For example, $\phi_{22} \approx 0$, but $\rho_X(2)$ is still large. This could suggest that $X_t$ only depends on $X_{t-1}$. In practice, the interpretation of sample ACF and PACF is not always as clear-cut.

## Autoregressive model of order 1

The AR(1) model is:

$$
X_t = \delta + \phi X_{t-1} + Z_t
$$

where $|\phi| <1$ (for stationarity) is an autoregressive coefficient, and $Z_t \overset{i.i.d.}{\sim} N(0, \sigma^2$ is WN which is independent of $X_{t-1}, X_{t-2}, \cdots, X_1$.

### Mean

The mean is given by

$$
\begin{gathered}
    E(X_t) = E(\sigma + \phi X_{t-1} + Z_t) = \sigma + \phi E(X_{t-1}) + E(Z_t) \\\\
    \mu = \sigma + \phi\mu \\\\
    E(X_t) = \mu = \frac{\sigma}{1 - \phi}
\end{gathered}
$$

From this we can easily see that when $\sigma = 0$, the mean is also 0.

An alternative expression for AR(1) is

$$
X_t - \mu = \phi(X_{t-1} - \mu) + Z_t
$$

This is the format R uses. The expectation is $E(X_t) = \mu$ because

$$
X_t = \underbrace{\mu - \phi\mu}\_{\sigma} + \phi X_{t-1} + Z_t
$$

### Variance

For variance and the ACF, we may set $\sigma = 0$ without losing any generalizability.

$$
\begin{aligned}
    Var(X_t) &= Var(\phi X_{t-1} + Z_t) \\\\
    &=Var(\phi X_{t-1}) + Var(Z_t) + \underbrace{2Cov(\phi X_{t-1}, Z_t)}\_{=0} \\\\
    &= \phi^2 Var(X_{t-1}) + \sigma^2
\end{aligned}
$$

Since $X_t$ is stationary, $Var(X_{t-1}) = Var(X_t)$. So,

$$
Var(X_t) = \frac{\sigma^2}{1 - \phi^2} = \gamma_X(0)
$$

We can also validate that the variance is always non-negative because $|\phi| < 1$. In R, we can use the `arima.sim` function to simulate an ARMA process.

### Autocovariance and autocorrelation

$$
\begin{aligned}
    \gamma_X(1) &= Cov(X_t, X_{t-1}) = Cov(\phi X_{t-1} + Z_t, X_{t-1}) \\\\
    &= \phi Cov(X_{t-1}, X_{t-1}) + Cov(Z_t, X_{t-1}) \\\\
    &= \phi Var(X_{t-1}) \\\\
    &= \phi \gamma_X(0) = \phi \frac{\sigma^2}{1 - \phi^2}
\end{aligned}
$$

From here we can find the autocorrelation with lag 1:

$$
\rho_X(1) = \frac{\gamma_X(1)}{\gamma_X(0)} = \phi
$$

Recall that in the figure, the first lag in the ACF plot is about 0.7.

Next we find the autocovariance with lag 2:

$$
\begin{aligned}
    \gamma_X(2) &= Cov(X_t, X_{t-2}) \\\\
    &= Cov(\phi X_{t-1} + Z_t, X_{t-2}) \\\\
    &= \phi Cov(X_{t-1}, X_{t-2}) + Cov(Z_t, X_{t-2}) \\\\
    &= \phi \gamma_X(1) = \phi^2 \frac{\sigma^2}{1 - \phi^2}
\end{aligned}
$$

Now $\rho_X(2)$ can be found:

$$
\rho_X(2) = \frac{\gamma_X(2)}{\gamma_X(0)} = \phi^2
$$

To generalize this, for $k \geq 1$, the autocovariance function is

$$
\begin{aligned}
    \gamma_X(k) &= Cov(X_t, X_{t-k}) \\\\
    &= Cov(\phi X_{t-1} + Z_t, X_{t-k}) \\\\
    &= \phi Cov(X_{t-1}, X_{t-k}) + Cov(\underbrace{Z_t, X_{t-k}}_{k \geq 1 \Rightarrow \text{indep.}}) \\\\
    &= \phi \gamma_X(k-1) \\\\
    &= \phi^2 \gamma_X(k-2) \\\\
    &= \cdots = \phi^k \gamma_X(0) = \phi^k \cdot \frac{\sigma^2}{1 - \phi^2}
\end{aligned}
$$

The autocorrelation function is

$$
\rho_X(k) = \frac{\gamma_X(k)}{\gamma_X(0)} = \phi^k
$$

We can see that $\rho_X(k)$ is an exponential function, and since $|\phi| < 1$, it exponentially decreases to 0 as the lag increases. In R, we can use the `ARMAacf` function to compute the theoretical ACF for an ARMA process. The same function can be used to find the PACF.

### Partial autocorrelation

Theoretically, The PACF is

$$
\phi_{00} = 1,\quad \phi_{11} = \phi = \rho(1), \quad \phi_{kk} = 0 \text{ for } k > 1
$$

The PACF of AR(1) cuts off after lag 1.
