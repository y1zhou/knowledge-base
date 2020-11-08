---
title: "Decomposition and Smoothing Methods"
date: 2020-11-02T11:12:20-05:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: ["Time Series", "Statistics", "Visualization", "R"] # list of related tags

slug: "time-series-decomposition-and-smoothing"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 60 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Smooth waves." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "2kWpWUHOheg" # Unsplash ID of the picture
---

The [reference text](https://online.stat.psu.edu/stat510/lesson/5/5.1) for this chapter is from Penn State University.

## Decomposition

In the previous chapters, differencing would be applied whenever we observe a mean trend or seasonal pattern. In this chapter we will learn another method -- `decomposition` procedures -- to extract the trend and seasonal factors in a time series using nonparametric detrending. Some more extensive decompositions might include long-run cycles, holiday effects, day of week effects etc.

One of the main objectives for a decomposition is to estimate seasonal effects that can be used to create and present `seasonally adjusted values`. A seasonally adjusted value removes the seasonal effect from a value so that trends can be seen more clearly.

### Basic structures

We're going to consider the following two structures for decomposition models:

1. Additive model:

    $$
    X_t = T_t + S_t + I_t = \text{Trend} + \text{Seasonal} + \text{Random}
    $$

2. Multiplicative model:

    $$
    X_t = T_t S_t I_t = \text{Trend} \times \text{Seasonal} \times \text{Random}
    $$

Here the random component is denoted $I_t$ because it's often called `irregular` in software.

The additive model is useful when the seasonal variation is relatively constant over time, whereas the multiplicative model is useful when the seasonal variation changes (increases) over time. This is because if we log-transform both sides of the multiplicative model, it becomes an additive model with a log-transformed series:

$$
\log X_t = \log T_t + \log S_t + \log I_t
$$

So when the variance is not constant over time, we can either apply a multiplicative model, or transform the series first and then apply an additive model. The latter approach is more flexible as it's not limited to a log-transformation.

### Steps in decomposition

1. Estimation of the trend. Two approaches with many variations of each are:

    - Smoothing techniques such as `moving averages`
    - Regression

2. Removal of the trend. The de-trending step is different for additive and multiplicative models:

    - Additive: $X_t - \hat{T}_t$
    - Multiplicative: $X_t / \hat{T}_t$

3. Estimation of seasonal components from the de-trended series. {{<hl>}}The seasonal effects are usually adjusted so that they average to 0 for an additive decomposition, or average to 1 for a multiplicative decomposition.{{</hl>}}

    The simplest method for estimating these effects is to average the de-trended values for a specific season. For example, we average the de-trended values for all Januarys in the series to get the seasonal effect for January.

4. Determination of the random (irregular) component.

    - Additive: $\hat{I}_t = X_t - \hat{T}_t - \hat{S}_t$
    - Multiplicative: $\hat{I}_t = \frac{X_t}{\hat{T}_t \times \hat{S}_t}$

    The remaining $\hat{I}_t$ part might be modeled with an ARIMA model. Sometimes steps 1 to 3 are iterated multiple times to get a stationary process. For example, after step 3 we could return to step 1 and estimate the trend based on the de-seasonalized series.

Let's say our times series can be decomposed into

$$
X_t = f(t) + W_t
$$

where $f(t)$ is a deterministic trend and $W_t$ is stationary. Our goal is to remove $f(t)$ from $X_t$.

If we take a parametric approach, we may assume some form of the deterministic trend such as

$$
f(t) = \beta_0 + \beta_1 t \text{ or } \beta_0 + \beta_1 t + \beta_2 t^2
$$

and the only thing we need to do is to estimate the $\beta$'s.

In our nonparametric approach[^nonparam], we estimate $f(t)$ directly at $t$. This is obviously more flexible as we don't have to make assumptions on the form of $f(t)$, but it comes with the caveat of being less efficient as we have to estimate $f(t)$ for each point of $t$. We also can't write out the formula of $f(t)$ explicitly -- we get $\hat{f}(t)$ instead.

[^nonparam]: The "nonparametric" here doesn't mean distribution-free. It's referring to the literal meaning of not having parameters.

Some popular nonparametric methods are kernel smoothing (LOWESS), splines, and wavelets. The kernel usually requires a larger sample size and the target function to be smooth. Splines are better with smaller sample sizes in higher dimensions, and wavelets are better with data that are discontinuous or contain spikes.

### Decomposition in R

The `decompose()` function in R can be used for both additive and multiplicative models. The `type` parameter can be specified as `"additive"` or `"multiplicative"`. Note that the seasonal span for a series needs to be defined using the `freq` parameter in the `ts()` function.

The following example uses the dataset of quarterly beer production in Australia from 1956 to 1973. Data for the first 18 years are used to demonstrate the functions. The series looks like this after additive decomposition[^beer-prod-r]:

{{< figure src="ausbeer_decomp_add.png" caption="Seasonally decomposed Australian beer production." numbered="true" >}}

[^beer-prod-r]: R code for generating the figure:

    ```r
    library(tidyverse)
    library(forecast)

    theme_set(
      ggpubr::theme_pubr() +
        theme(
          axis.title = element_text(face = "bold", size = 14)
        )
    )

    data("ausbeer", package = "fpp")
    beerprod <- window(ausbeer, end = c(1973, 4))
    decomp_beer_add <- decompose(beerprod, type = "additive")
    autoplot(decomp_beer_add, range.bars = F)
    ```

The function plots the original data, the mean trend, the seasonal component, and the remainders[^remainder-vs-random]. See that the seasonal component stays constant for each quarter and is a regularly repeating pattern. We can also set the `type` parameter to "multiplicative" and get:

[^remainder-vs-random]: The word "remainder" is used rather than "random" because maybe it's not really random!

{{< figure src="ausbeer_decomp_multi.png" caption="Australian beer production series after multiplicative decomposition." numbered="true" >}}

### Smoothing

Smoothing is used to help us better see patterns in time series. Generally we try to smooth out the irregular roughness to get a cleaner signal. For seasonal data, we might also smooth out the seasonality so that the trend could be identified. Smoothing doesn't provide us with a model, but it can be a good step in describing various components of a series.

The term `filter` is sometimes used to describe a smoothing procedure. For example, if the smoothed value for a particular time is calculated as a linear combination of the observations for surrounding time points, it might be said that we've applied a `linear filter` to the data.

#### Moving average

As we said earlier, a simple but effective method to estimate the trend is the moving average. For example, the trend values in the Australian beer example were determined as centered moving averages with length 4 (because there are four quarters per year). For time $t=3$, we first average the observed data values at times 1 to 4:

$$
\frac{1}{4} (X_1 + X_2 + X_3 + X_4)
$$

Then the average of the values at times 2 to 5:

$$
\frac{1}{4}(X_2 + X_3 + X_4 + X_5)
$$

Then we take the average of these two values to get

$$
\frac{1}{8}X_1 + \frac{1}{4}X_2 + \frac{1}{4}X_3 + \frac{1}{4}X_4 + \frac{1}{8}X_5
$$

We're doing this because the length is an even number. For length 5, we can simply take

$$
\frac{1}{5}(X_1 + X_2 + X_3 + X_4 + X_5)
$$

More generally, the `centered moving average` for time $t$ of a quarterly series is

$$
\frac{1}{8}X_{t-2} + \frac{1}{4}X_{t-1} + \frac{1}{4}X_t + \frac{1}{4}X_{t+1} + \frac{1}{8}X_{t+2}
$$

and for a monthly series:

$$
\frac{1}{24}X_{t-6} + \sum_{j=-5}^5 \frac{1}{12}X_{t+j} + \frac{1}{24}X_{t+6}
$$

The `decompose()` function uses this method (symmetric window with equal weight), and we can get the results with `decomp_beer_add$trend`. See how we have two `NA` values at the beginning and at the end because the moving average can't be calculated at these time points. This linear filter can be calculated manually in R using:

```r
trend_pattern <- filter(beerprod, sides = 2,
                        filter = c(1/8, 1/4, 1/4, 1/4, 1/8))
```

Other variations of the moving average exist and may serve different purposes. For example, to remove a seasonal effect, we may calculate the moving average with a length equal to the seasonal span. We may only use the previous $n$ time points, i.e. the one-sided moving average:

$$
\frac{1}{4}(X_t + X_{t-1} + X_{t-2} + X_{t-3})
$$

In R, this can be found by

```r
filter(beerprod, filter = rep(1/4, 4), sides = 1)
```

#### Lowess

Another more advanced decomposition method is the locally weighted scatterplot smoothing (LOWESS)[^lowess]. To understand what this is, we need to introduce the concept of kernel smoothing.

[^lowess]: For a more in-depth introduction to Lowess, see [this chapter in the Nonparametric Methods course]({{< relref "/series/nonparam-stat/5-modern-methods/5.3-nonparametric-regression/index.md" >}}).

In the simplest univariate case where $X_1, \cdots, X_n \overset{i.i.d.}{\sim} f_X$, our goal is to estimate the density function $f_X$. In parametric density estimation, we assume some density function (e.g. Binomial, Poisson, Normal, etc.) and estimate the parameters, often with MLE.

Parametric density estimation is efficient when the assumptions on the distribution are met. In nonparametric density estimation, we estimate $\hat{f}_X(x)$ at a given $x$ without assuming a form of the distribution. At the expense of efficiency, we gain a lot of flexibility.

The histogram is a naive approach of density estimation, although an obvious disadvantage is that it's discrete and not smooth. It also depends on the bin width $b$ and the starting point. `Kernel smoothing` is proposed to overcome some of these problems.

Let $h$ be the bandwidth. We can express the density estimate at $x$ of the histogram as

$$
\begin{aligned}
    \hat{f}_X(x) &= \frac{\text{Frequency}}{nb} \\\\
    &= \frac{\text{number of }X_i \in (x-h, x+h)}{2nh} \\\\
    &= \frac{\sum\_{i=1}^n I(X_i \in (x-h, x+h))}{2nh} \\\\
    &= \frac{1}{2nh}\sum\_{i=1}^n I\left( |x - X_i| < h \right) \\\\
    &= \frac{1}{2nh}\sum\_{i=1}^n I\left( \left|\frac{x - X_i}{h}\right| < 1 \right)
\end{aligned}
$$

The discreteness comes from the indicator function. The idea of kernel smoothing is instead of this indicator function, we plug in a smoother `kernel function` $K$ that assigns larger weights for points closer to $x$ and lower weights for distant points:

$$
\hat{f}(x) = \frac{1}{2nh}\sum_{i=1}^n K\left( \frac{x - X_i}{h} \right)
$$
