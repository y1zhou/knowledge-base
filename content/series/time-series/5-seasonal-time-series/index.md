---
title: "Seasonal Time Series"
date: 2020-10-13T23:36:46-04:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: ["Time Series", "R"] # list of related tags

slug: "seasonal-time-series"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 50 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Aerial view of a forest road in autumn." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "Qy-CBKUg_X8" # Unsplash ID of the picture
---

In the last chapter, we learned about how to deal with mean and variance changes in time series data. Periodic patterns are also frequently observed, and we will use seasonal ARIMA models for this type of data.

## Motivating example

The data are the monthly sales of six types of wines Australian wines (red, rose, sweet white, dry white, sparkling, and fortified) for years 1980-1994. The units are thousands of liters. We'll focus on the red wine sales data, and use the `ts()` function to turn the column into a time series object.

```r
wine <- data.table::fread("https://raw.githubusercontent.com/y1zhou/knowledge-base/master/content/series/time-series/5-seasonal-time-series/wine.csv")
wine <- ts(data = wine$Red, start = 1980, end = c(1994, 12), frequency = 12)
```

As always, we first plot the series and its ACF and PACF. On top of the mean and variance increases over time, we also observe a periodical pattern that repeats every year. This means we might be able to use the sales in June 1985 to predict the sales in June 1986. In seasonal AR(1) (MA(1)), $X_{t-12}$ ($Z_{t-12}$) can be used to predict $X_t$; in seasonal AR(2), $X_{t-12}$ and $X_{t-24}$ can be used to predict $X_t$, etc.

{{< figure src="wine_ts.png" caption="Time series plot, ACF and PACF of the wine dataset." numbered="true" >}}

When we see a seasonal pattern in the time series plot, it's often desirable to increase the `lag.max` parameter when plotting the ACF and PACF to see the entire pattern. The ACF slowly decreases as a result of the series being trended. Because the series is seasonal, the ACF for the seasonal lags (12, 24, 36, etc. multiples of the seasonal frequency) are higher than the ACF of the other lags. We see a combination of these effects for this trended, seasonal data.

This effect could be removed using `seasonal differencing`. But before that, let's see some simpler models first.

## Seasonal ARMA

A seasonal pattern occurs when a time series is affected by seasonal factors such as the time of the year or the day of the week. {{<hl>}}Seasonality is always of a fixed and known period{{</hl>}}. The `frequency` is the number of observations before the seasonal pattern repeats[^frequency-def].

[^frequency-def]: This is the opposite of the definition of frequency in physics, where this would be called the period and its inverse would be called the frequency.

### Seasonal MA(1)

We will start with the `seasonal MA(1)` model. As we have monthly data in the motivating example, let's say the `period` $s=12$. The $MA(1)_{s=12}$ model is

$$
X_t = \mu + Z_t + \Theta Z_{t-12}, \quad Z_t \overset{i.i.d.}{\sim} WN(0, \sigma^2)
$$

Recall that the MA(1) model had a linear combination of $Z_t$ and $Z_{t-1}$, but in our seasonal MA(1) model we have a combination of $Z_t$ and $Z_{t-12}$ where the 12 corresponds to the period. If we didn't know there's a seasonal pattern, we could call this an MA(12) model. The $\Theta$ notation is to emphasize the seasonality.

Calculating the ACF is similar to the process for the MA(1) model. We first find the autocovariance function:

$$
\begin{aligned}
    \gamma_X(0) &= Var(X_t) = Var(Z_t) + \Theta^2 Var(Z_{t-12}) + 0 \\\\
    &= (1 + \Theta^2)\sigma^2 \\\\
    \gamma_X(k) &= Cov(X_t, X_{t-k}), \quad k \geq 1 \\\\
    &= Cov(Z_t + \Theta Z_{t-12}, X_{t-k}) \\\\
    &= \underbrace{Cov(Z_t, X_{t-k})}_{=0} + \Theta Cov(Z_{t-12}, X_{t-k}) \\\\
    &= \Theta Cov(Z_{t-12}, Z_{t-k} + \Theta Z_{t-12-k}) \\\\
    &= \Theta Cov(Z_{t-12}, Z_{t-k}) + \Theta^2 \underbrace{Cov(Z_{t-12}, Z_{t-12-k})}_{=0} \\\\
    &= \begin{cases}
        \Theta \sigma^2, & k = 12 \\\\
        0, & \text{otherwise}
    \end{cases}
\end{aligned}
$$

Thus the ACF is

$$
\rho_X(k) = \begin{cases}
    \frac{\Theta}{1 + \Theta^2}, & k = 12 \\\\
    0, & \text{otherwise}
\end{cases}
$$

To find the theoretical ACF, we may use the `ARMAacf()` function in R. The theoretical PACF can be found by setting parameter `pacf = TRUE`. We can see that all ACF values are zeros except at $k=12$. For the PACF, there are spikes at lags that are multiples of 12.

```r
ARMAacf(ma = c(rep(0, 11), -0.5), lag.max = 30)
#    0    1    2    3    4    5    6    7    8    9   10
#  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#   11   12   13   14   15   16   17   18   19   20   21
#  0.0 -0.4  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
#   22   23   24   25   26   27   28   29   30
#  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

ARMAacf(ma = c(rep(0, 11), -0.5), lag.max = 50, pacf = T)
#  [1]  0.00000000  0.00000000  0.00000000  0.00000000
#  [5]  0.00000000  0.00000000  0.00000000  0.00000000
#  [9]  0.00000000  0.00000000  0.00000000 -0.40000000
# [13]  0.00000000  0.00000000  0.00000000  0.00000000
# [17]  0.00000000  0.00000000  0.00000000  0.00000000
# [21]  0.00000000  0.00000000  0.00000000 -0.19047619
# [25]  0.00000000  0.00000000  0.00000000  0.00000000
# [29]  0.00000000  0.00000000  0.00000000  0.00000000
# [33]  0.00000000  0.00000000  0.00000000 -0.09411765
# [37]  0.00000000  0.00000000  0.00000000  0.00000000
# [41]  0.00000000  0.00000000  0.00000000  0.00000000
# [45]  0.00000000  0.00000000  0.00000000 -0.04692082
# [49]  0.00000000  0.00000000
```

### Seasonal MA(Q)

In general the $MA(Q)_s$ model is

$$
X_t = Z_t + \Theta_1 Z_{t-s} + \Theta_2 Z_{t-2s} + \cdots + \Theta_Q Z_{t-Qs}, \quad Z_t \overset{i.i.d.}{\sim} WN(0, \sigma^2)
$$

The ACF has spikes at lags $s, 2s, \cdots, Qs$ and is zero afterwards. The PACF exponentially decreases at lags that are multiples of $s$ and is zero for other lags.

### Seasonal AR(1)

The $AR(1)_{12}$ model is

$$
X_t = \Phi X_{t-12} + Z_t, \quad Z_t \sim WN(0, \sigma^2)
$$

Again, this is similar to the AR(1) model with the $X_{t-1}$ term replaced by $X_{t-12}$. The ACF is

$$
\rho_X(12k) = \Phi^k, \quad k = 0, 1, 2, \cdots
$$

which is found from

$$
\begin{aligned}
    \gamma_X(0) &= \frac{\sigma^2}{1 - \Phi^2} \\\\
    \gamma_X(k) &= \begin{cases}
        \frac{\sigma^2\Phi^k}{1 - \Phi^2}, & k = 12, 24, 36, \cdots \\\\
        0, & \text{otherwise}
    \end{cases}
\end{aligned}
$$

### Seasonal AR(P)

The $AR(P)_s$ model is

$$
X_t = \Phi_1 X_{t-s} + \Phi_2 X_{t-2s} + \cdots + \Phi_P X_{t-Ps} + Z_t, \quad Z_t \sim WN(0, \sigma^2)
$$

The ACF exponentially decreases at lags that are multiples of $s$. The PACF spikes at lags that are multiples of $s$ and takes value zero for all other lags.

### Seasonal ARMA(1, 1)

The $ARMA(1, 1)_{12}$ model is

$$
X_t = \Phi X_{t-12} + Z_t - \Theta Z_{t-12}, \quad Z_t \sim WN(0, \sigma^2)
$$

Alternatively, we can write it using the backshift operator:

$$
(1 - \Phi B^{12})X_t = (1 - \Theta B^{12})Z_t
$$

The conditions for **stationarity and invertibility** are the same as the ones in the regular ARMA model: the absolute values of the roots of the characteristic equations should be larger than 1.

$$
\begin{gathered}
    \text{Stationarity: } |\Phi| > 1 \\\\
    \text{Invertibility: } |\Theta| > 1
\end{gathered}
$$

## Multiplicative seasonal ARMA

Now we are going to combine the regular ARMA and seasonal ARMA models. This is done by multiplying the AR/MA polynomials. For example, the multiplicative combination of $MA(1)$ and $MA(1)\_{12}$ is denoted $ARMA(0, 1) \times (0, 1)\_{s=12}$, and can be expressed as

$$
X_t = (1 + \theta B)(1 + \Theta B^{12})Z_t = Z_t + \theta Z_{t-1} + \Theta Z_{t-12} + \theta\Theta Z_{t-13}
$$

It's not hard to show that $\gamma(0) = (1 + \theta^2)(1 + \Theta^2)\sigma^2$. We can also show that the ACF would have spikes at not only lags 1 and 12, but also 11 and 13.

$$
\begin{gathered}
    \rho(1) = \frac{\theta}{1 + \theta^2}, \quad \rho(12) = \frac{\Theta}{1 + \Theta^2}, \\\\
    \rho(11) = \rho(13) = \frac{\theta\Theta}{(1 + \theta^2)(1 + \Theta^2)}
\end{gathered}
$$

The general form of an $ARMA(p, q) \times (P, Q)_s$ model is

$$
\Phi(B)\phi(B)X_t = \Theta(B)\theta(B)Z_t, \quad Z_t \overset{i.i.d.}{\sim}WN(0, \sigma^2)
$$

where

$$
\begin{gathered}
    \phi(x) = 1 - \phi_1 x - \phi_x x^2 - \cdots - \phi_p x^p, \\\\
    \Phi(x) = 1 - \Phi_1 x^s - \Phi_x x^{2s} - \cdots - \Phi_P x^{Ps}, \\\\
    \theta(x) = 1 + \theta_1 x + \theta_2 x^2 + \cdots + \theta_q x^q, \\\\
    \Theta(x) = 1 + \Theta_1 x^s + \Theta_2 x^{2s} + \cdots + \Theta_Q x^{Qs} \\\\
\end{gathered}
$$

Here the AR order is $p + Ps$, the MA order is $q + Qs$, and the number of parameters is $p + P + q + Q$. We can also include an intercept term and add one to the total number of parameters. **The order of the seasonal part of the model rarely goes higher than 1**, especially after we include differencing.

For example, $ARMA(1, 0) \times (0, 2)_{12}$ can be written as

$$
(1 - \phi B) X_t = (1 + \Theta_1 B^{12} + \Theta_2 B^{24}) Z_t,
$$

and the following model

$$
X_t = 0.6 X_{t-12} + Z_t + 0.4 Z_{t-1} - 0.8 Z_{t-12} - 0.32 Z_{t-13}
$$

is $ARMA(0, 1) \times (1, 1)_{12}$.

### ACF and PACF
