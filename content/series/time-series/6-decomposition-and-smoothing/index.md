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

The following example uses the dataset of quarterly beer production in Australia from 1956 to 1973. Data for the first 18 years are used to demonstrate the functions. The decomposed series looks like this[^beer-prod-r]:

{{< figure src="ausbeer_decomp.png" caption="Seasonally decomposed Australian beer production." numbered="true" >}}

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
    decomp_beer <- decompose(beerprod, type = "additive")
    autoplot(decomp_beer, range.bars = F)
    ```
