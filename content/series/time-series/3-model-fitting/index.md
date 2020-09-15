---
title: "Model Fitting and Forecasting"
date: 2020-09-14T12:06:41-04:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: ["Statistics", "Time Series", "Regression", "R"] # list of related tags

slug: "time-series-model-fitting-and-forecasting"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 30 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

In the previous chapter, we learned about the characteristics of AR, MA and ARMA processes and generated data from these known models. In real data analysis, we only have the data and don't know the true underlying model.

This model-building strategy consists of three steps: model specification (identification), model fitting, and model diagnostics.

## Model identification

In model identification, we need to find appropriate models for a given observed series. The model is tentative, and is subject to revision later on in the analysis.

Usually the first step is plotting the data and see if it looks stationary. If it seems stationary, we try to identify serial correlation using the ACF and the PACF. The characteristics of the ACF and PACF are given in the table below.

| Process      | ACF                            | PACF                           |
| ------------ | ------------------------------ | ------------------------------ |
| $AR(p)$      | Tails off as exponential decay | Cuts off after lag $p$         |
| $MA(q)$      | Cuts off after lag $q$         | Tails off as exponential decay |
| $ARMA(p, q)$ | Tails off after lag $q-p$      | Tails off after lag $p-q$      |

### Information criteria

For each model, we start from the lowest-order model[^model-identification-principle] and look at the residuals. If we observe prominent spikes in either the ACF or the PACF, we then sequentially increase the order. Otherwise, we keep the model as simple as possible.

[^model-identification-principle]: Follow the [principle of parsimony](https://en.wikipedia.org/wiki/Occam%27s_razor).

Two popular criteria to choose the best model are the `AIC` (Akaike information criterion) and the `BIC` (Bayesian information criterion):

$$
\begin{gathered}
    AIC(\ell) = \ln(\tilde{\sigma}\_\ell^2) + \frac{2(\ell + 1)}{T} \\\\
    BIC(\ell) = \ln(\tilde{\sigma}\_\ell^2) + \frac{(\ell + 1)\ln(T)}{T}
\end{gathered}
$$

Both measures were mentioned in regression analysis. They deal with the trade-off between the goodness of fit of the model and the simplicity of the model, i.e. risks of overfitting and underfitting.

Recall the definition of the [likelihood function]({{< relref "/series/maths-stat/9-estimation-under-parametric-models/9.1-maximum-likelihood-estimator/index.md#maximum-likelihood-estimator" >}}). It's a measure of how likely the parameter estimates are.

## Model fitting

We may fit 2-3 suitable models and check their assumptions...
