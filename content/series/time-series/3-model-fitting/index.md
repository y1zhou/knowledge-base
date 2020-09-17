---
title: "Model Fitting and Forecasting"
date: 2020-09-14T12:06:41-04:00
summary: "This model-building strategy consists of three steps: model specification (identification), model fitting, and model diagnostics." # appears in list of posts
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

In the previous chapter, we learned about the characteristics of AR, MA and ARMA processes and generated data from these known models. In real data analysis, we only have the data and don't know the true underlying model. Our workflow typically follows:

1. Time series plot: stationarity (mean, variance trend, seasonal pattern).
2. Identification of serial correlation using ACF and PACF.
3. Model fitting.
4. Model diagnosis.
5. Choose the best model by the AIC/BIC (minimize), prediction ability (MAPE), and the number of parameters.

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

Recall the definition of the [likelihood function]({{< relref "/series/maths-stat/9-estimation-under-parametric-models/9.1-maximum-likelihood-estimator/index.md#maximum-likelihood-estimator" >}}). It's a measure of how likely the parameter estimates are. The problem with likelihood functions is that they always increase as we use more parameters. We introduce information criteria to deal with the trade-off between the goodness of fit of the model and the simplicity of the model, i.e. risks of overfitting and underfitting.

Two popular criteria to choose the best model are the `AIC` (Akaike information criterion) and the `BIC` (Bayesian information criterion):

$$
\begin{gathered}
    AIC(l) = \ln(\tilde{\sigma}\_l^2) + \frac{2(l + 1)}{T} \\\\
    BIC(l) = \ln(\tilde{\sigma}\_l^2) + \frac{(l + 1)\ln(T)}{T}
\end{gathered}
$$

where $\tilde{\sigma}\_l^2$ is the maximum likelihood estimate of $\sigma^2$, $T$ is the sample size, and $l$ is the number of parameters[^num-of-parameters]. We select $l \in \\{0, 1, \cdots, P\\}$ that minimizes the AIC or the BIC. Essentially we're trying to minimize the negative log-likelihood considering the number of parameters.

[^num-of-parameters]: $l+1$ puts an emphasize on the intercept term in the model.

To determine whether or not to include the intercept in the model, we may use a one-sample t-test with hypotheses $H_0: \mu = 0$ vs. $H_1: \mu \neq 0$, and the test statistic is

$$
t = \frac{\bar{X}}{\sqrt{\frac{Var(X_t)}{T}}} \approx N(0, 1)
$$

Here we're essentially investigating $\mu = E(X_t)$.

## AR

In this section, we're going to use the AR(1) process as an example to introduce parameter estimation and model dignostics methods.

### Estimation

The common methods to estimate parameters are `MME` (method of moments estimation), `LSE` (least squares estimation) and `MLE` (maximum likelihood estimation), where the MLE has a distribution assumption as the likelihood function requires the PDF. In many cases the three methods are very similar, and the MLE performs better in complex cases.

#### Method of moments

For an AR(1) process with mean 0, the model is given by

$$
X_t = \phi X_{t-1} + Z_t
$$

where $Z_t: WN(0, \sigma^2)$. The parameters to estimate are $\phi$ and $\sigma^2$.

$$
\begin{gathered}
    \phi = \frac{\gamma_X(1)}{\gamma_X(0)}, && \sigma^2 = (1 - \phi^2)\gamma_X(0) \\\\
    \hat\phi = \frac{\sum_{t=1}^{T-1} X_t X_{t+1}}{\sum_{t=1}^T X_t^2}, && \hat\sigma^2 = (1 - \hat\phi^2) \cdot \frac{1}{T}\sum_{t=1}^T X_t^2
\end{gathered}
$$

For an AR(1) model with mean $\mu$, the model is

$$
X_t - \mu = \phi(X_{t-1} = \mu) + Z_t
$$

so we have one more parameter to estimate. The formulae are very similar:

$$
\begin{gathered}
    \hat\mu = \bar{X} \\\\
    \hat\phi = \frac{\sum_{t=1}^{T-1} (X_t - \bar{X})(X_{t+1} - \bar{X})}{\sum_{t=1}^T (X_t - \bar{X})^2} \\\\
    \hat\sigma^2 = (1 - \hat\phi^2) \cdot \frac{1}{T}\sum_{t=1}^T (X_t - \bar{X})^2
\end{gathered}
$$

For an AR(2) model with mean 0:

$$
X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + Z_t,
$$

we may use the Yule-Walker equation to get

$$
\begin{gathered}
    \rho_X(1) = \phi_1 + \phi_2 \rho_X(1) \\\\
    \rho_X(2) = \phi_1 \rho_X(1) + \phi_2
\end{gathered}
$$

which allows us to express $\phi_1$ and $\phi_2$ using the autocorrelation functions:

$$
\begin{gathered}
    \phi_2 = \frac{\rho_X(2) - \rho_X(1)^2}{1 - \rho_X(1)^2} \\\\
    \phi_1 = \frac{\rho_X(1) - \rho_X(1) \rho_X(2)}{1 - \rho_X(1)^2}
\end{gathered}
$$

We know that

$$
\hat\rho_X(k) = \frac{\sum_{t=1}^{T_k} X_t X_{t+k}}{\sum_{t=1}^T X_t^2},
$$

therefore, $\hat\phi_1$ and $\hat\phi_2$ can both be expressed using $\hat\rho_X(1)$ and $\hat\rho_X(2)$. The estimator for $\sigma^2$ can be found from the equation for $Var(X_t)$:

$$
\hat\sigma^2 = \left(1 - \hat\phi_1 \hat\rho_X(1) - \hat\phi_2 \hat\rho_X(2) \right) \cdot \frac{1}{T}\sum_{t=1}^T X_t^2
$$

#### Least squares

For the AR(1) process

$$
X_t - \mu = \phi(X_{t-1} - \mu) + Z_t,
$$

the sum of sqaured errors is

$$
SSE = \sum_{t=1}^T Z_t^2 = \sum_{t=2}^T (X_t - \delta - \phi X_{t-1})^2
$$

To minimize the conditional sum of squares:

$$
\begin{gathered}
    s^*(\phi, \delta) = \sum_{t=2}^T \left( (X_t - \mu) - \phi(X_{t-1} - \mu) \right)^2 \\\\
    \hat\mu \approx \bar{X} \\\\
    \hat\phi = \frac{\sum_{t=2}^T (X_t - \bar{X})(X_{t-1} - \bar{X})}{\sum_{t=2}^T (X_{t-1} - \bar{X})^2} \\\\
    \hat\sigma^2 = \frac{1}{T-3}\sum_{t=2}^T \left( (X_t - -\bar{X}) - \hat\phi(X_{t-1} - \bar{X}) \right)^2
\end{gathered}
$$

#### Maximum likelihood estimation

Assuming $Z_t \overset{i.i.d.}{\sim} N(0, \sigma^2)$, we want to find $(\mu, \phi, \sigma^2)$ which minimizes

$$
L(\phi, \mu, \sigma^2) - \prod_{t=2}^T \frac{1}{\sqrt{2\pi}\sigma}\exp \left\\{ -\frac{((X_t - \mu) - \phi(X_{t-1} - \mu))^2}{2\sigma^2} \right\\}
$$

By partial derivatives we'll get

$$
\begin{gathered}
    \hat\mu = \frac{1}{T-1}\sum_{t=2}^T X_{t-1} \approx \bar{X} \\\\
    \hat\phi = \frac{\sum_{t=2}^T (X_t - \hat\mu)(X_{t-1} - \hat\mu)}{\sum_{t=2}^T (X_{t-1} - \hat\mu)^2} \\\\
    \hat\sigma^2 = \frac{1}{T-1}\sum_{t=2}^T \left( (X_t - \hat\mu) - \hat\phi(X_{t-1} - \hat\mu) \right)^2
\end{gathered}
$$

### Model diagnostics

To check the assumption that $Z_t \sim N(0, \sigma^2)$, we perform the `white noise check`.

The standardized residuals (innovations) are

$$
e_t = \frac{X_t - \tilde{X}_t^{t-1}}{\sqrt{\hat{V} \left(X_t - \tilde{X}_t^{t-1} \right)}}
$$

where $\tilde{X}_t^{t-1}$ is the forecast for time $t$ using information up to time $t-1$, and $\hat{V}(X_t - \tilde{X}_t^{t-1})$ is the estimate of $Var(X_t - \tilde{X}_t^{t-1})$. If the correct model is chosen, these standardized residuals should approximately follow $N(0, 1)$. Thus, we'd expect about 95% of them to fall within $\pm 2$ and about 99% to fall within $\pm 3$. Standardized residulas falling outside of this range indicate an incorrect model may have been chosen.

We should also check the ACF/PACF of the standardized residuals. If the residuals are really i.i.d. normal, there should only be a spike at lag 0 and all the other lags should fall within the band:

$$
\hat\rho(k) \pm \frac{2}{\sqrt{T}}
$$

This could be tested using the `Ljung-Box test` (Portmanteau test), which tests if any of a group of autocorrelations of a time series are different from zero. The hypotheses are

$$
\begin{gathered}
    H_0: \rho(1) = \cdots = \rho(m) = 0 \text{ vs.} \\\\
    H_1: \rho(i) \neq 0 \text{ for at least one } i.
\end{gathered}
$$

The test statistic is

$$
Q(m) = T(T+2) \sum_{k+1}^m \frac{\hat\rho(k)^2}{T-k} \approx \chi^2(m - \text{ number of parameters in the model})
$$

Usually we have $m \approx \log(T)$ where $T$ is the number of timepoints. For AR(1) with mean 0, there's only one parameter[^box-test-df] $\phi$ so the degrees of freedom is $m-1$. For AR(1) with mean $\mu$, the d.f. is $m-2$.

[^box-test-df]: We ignore $\sigma^2$ when counting the number of parameters.

#### Simulation in R

We use the following code to simulate an AR(1) process with $\phi = 0.5$ and $n = 100$.

```r
set.seed(1)
x <- arima.sim(list(order = c(1, 0, 0), ar = 0.5), n = 100)
x.fit <- arima(x, order = c(1, 0, 0), include.mean = F)
x.fit
# Call:
#   arima(x = x, order = c(1, 0, 0), include.mean = F)
#
# Coefficients:
#   ar1
# 0.5377
# s.e.  0.0872
#
# sigma^2 estimated as 0.8398:  log likelihood = -133.33,  aic = 270.67
```

The estimated parameter $\hat\phi = 0.5377$ and the estimated $\hat\sigma^2 = 0.8398$. Thus the fitted model is

$$
X_t = 0.5377 X_{t-1} + Z_t, \quad Z_t \sim N(0, 0.8398)
$$

Next we check the assumptions of the fitted model. The `tsdiag()` function generates three diagnosis plots: the standardized residuals plot, the ACF of the residuals, and the p-values for the Ljung-Box statistic. The `gof.lag` parameter gives the maximum number of lags for the test.

Alternatively, we may use functions from other packages to model (`forecast::auto.arima`) and plot (`ggfortify::ggtsdiag`) the data[^ggfortify].

[^ggfortify]: The function returns a `ggmultiplot` class object.

    ```r
    library(ggfortify)
    ggtsdiag(x.fit, gof.lag = 20) +
      ggpubr::theme_pubr()
    ```

{{< figure src="ts_diag.png" caption="Diagnosis plots for the simulated AR(1) series." numbered="true" >}}

The standardized residuals plot is centered around zero, and there's no specific patterns with most values falling within $\pm 2$. For the ACF plot, everything falls within the bands except for lag 0, indicating no serial correlations. Finally, the p-values for the Ljung-Box statistic are all very high, meaning the none of the $H_0$'s can be rejected and there's no evidence at any lag that $\rho = 0$. Overall we observe no signs of violation of the assumptions.

We can also perform the Ljung-Box test using the `Box.test()` function. Here's an example for lag 3:

```r
Box.test(x.fit$residuals, lag = 3, type = "Ljung-Box", fitdf = 1)
# 	Box-Ljung test

# data:  x.fit$residuals
# X-squared = 1.5143, df = 2, p-value = 0.469
```

Note we specified `fitdf = 1` to subtract 1 from the degrees of freedom to get the correct number. The d.f. should be 2, which is $m$ (lag 3) minus the number of parameters in AR(1). This is why the computed p-value of 0.469 doesn't seem to match the value in the figure[^ljung-box-df].

[^ljung-box-df]: The degrees of freedom in R's function is wrong. It's fixed at $m$.