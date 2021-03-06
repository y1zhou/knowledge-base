---
title: "Penalized Linear Regression and Model Selection"
date: 2021-04-12T10:20:00-05:00
summary: 'This lecture covers some Bayesian connections to penalized regression methods such as ridge regression and the LASSO. Further discussion of the posterior predictive distribution as well as model selection criterion (DIC) is included.' # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-penalized-linear-regression-and-model-selection"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 110 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Lego heads." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "AVmmam25Mpg" # Unsplash ID of the picture
---



This lecture discusses some further considerations with regards to Bayesian linear regression. Penalized regression methods can be thought of as applying "shrinkage" to regression parameters towards zero to improve model fit and prediction.

Note that this is only meant to be an introduction of these topics. For a more comprehensive overview, *Bayesian Data Analysis* by Gelman et al. is a good reference.


## Diabetes example

The motivating example used comes from *Hoff*, page 161. Consider a sample of `\(n=442\)` diabetes patients for whom a measure of disease progression after one year, `\(y\)`, is recorded along with baseline measurements for ten other variables `\(x_1, \cdots, x_{10}\)`.

Since it is believed the relationship is nonlinear with numerous interactions, we consider a total of `\(p=64\)` predictor variables, obtained from:

- main effects of the ten variables,
- interactions of the form `\(x_i x_j\)` for all `\(i \neq j\)`, and
- quadratic terms `\(x_i^2\)`.

We will fit regression models using a training sample of `\(n_{\text{train}} = 342\)` patients. The model is then tested by evaluating how well it predicts disease progression on the remaining test sample of `\(n_{\text{test}}=100\)` patients.


```r
library(tidyverse)
library(ggpubr)
library(rjags)
library(MCMCvis)
```



```r
load("diab.Rdata")

X_train <- diab_train %>%
  select(-y) %>%
  as.matrix()

X_test <- diab_test %>%
  select(-y.te) %>%
  as.matrix()
```


## Bayesian linear regression

Suppose we are interested in regressing `\(Y_i\)` on predictor variables `\(x_{i1}, x_{i2}, \cdots, x_{ip}\)` for `\(i = 1, \cdots, n\)`. The **sampling model** in our `Bayesian linear regression model` looked like this:

<div>
$$
Y_i \mid \mu_i, \sigma^2 \overset{\text{indep.}}{\sim} \mathcal{N}(\mu_i, \sigma^2), \quad i = 1, \cdots, n
$$
</div>

where

`\begin{equation}
\mu_i = \beta_0 + \beta_1 x_{i1} + \cdots + \beta_p x_{ip} = \beta_0 + \sum_{j=1}^p \beta_j x_{ij} \tag{1}
\end{equation}`

Here the **prior models** are:

`\begin{gathered}
\beta_0, \cdots, \beta_p \overset{iid}{\sim} \mathcal{N}(0, \kappa^2) \text{ for some large } \kappa > 0 \\
\sigma^2 \sim \text{IG}(0.01, 0.01)
\end{gathered}`

Nothing new so far!


### High-dimensional linear regression

Some issues can arise in linear regression when we have a large number (`\(p\)`) of predictor variables.

First, many of the predictor variables may have minimal impact on explaining the response variable. For example, in Eq. (1) many of the `\(\beta_i\)`'s may be equal or close to zero.

Second, there is no **unique** optimal set of regression coefficients when the number of predictors `\(p\)` is larger than the sample size `\(n\)`. Recall that in the ordinary least squares formulation, the OLS regression coefficients estimate `\(\hat\beta_{\text{OLS}}\)` minimizes the following:

<div>
$$
\sum_{i=1}^n \left[ y_i - \left( \beta_0 + \sum_{j=1}^p \beta_j x_{ij} \right) \right]^2 = (\boldsymbol{y} - \boldsymbol{X\beta})^\prime (\boldsymbol{y} - \boldsymbol{X\beta})
$$
</div>

The solution in matrix notation is:

<div>
$$
\hat\beta_{\text{OLS}} = (\boldsymbol{X}^\prime \boldsymbol{X})^{-1} \boldsymbol{X}' \boldsymbol{y}
$$
</div>

When `\(p > n\)`, `\(\boldsymbol{X}'\boldsymbol{X}\)` is not full rank and is thus not invertible.

The third problem is `overfitting`. We may have a model built on an initial set of training data which fits very well, but does not *generalize* to new data well. This is known as the `Bias-variance trade-off`.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/bias-var-tradeoff-1.svg" alt="Demonstration of the bias-variance trade-off." width="384" />
<p class="caption">Figure 1: Demonstration of the bias-variance trade-off.</p>
</div>

In the illustration above, the green model is "simple", which means it's biased but also has less variability in the estimated line. The red model is complex and not biased at all, but it's highly variable depending on the observed training data. When given a new data set, we might trust the green model more than the red one as it generalizes better.


### Penalized linear regression

Penalized regression methods address these three issues by estimating **biased** regression coefficients in exchange for lower variance, and generally improves model fit and predictive performance. Common penalized regression methods include:

- Ridge regression
- LASSO
- Elastic net regression -- a combination of ridge and LASSO

All of these involve adding a **penalty term** to the function being minimized, i.e. the optimal regression coefficients `\(\hat\beta_{\text{pen}}\)` minimize:

<div>
$$
\sum_{i=1}^n \left[ y_i - \left( \beta_0 + \sum_{j=1}^p \beta_j x_{ij} \right) \right]^2 + \text{Pen}(\boldsymbol{\beta})
$$
</div>

This penalty term differs in the three methods.


#### Ridge regression

Ridge regression was introduced by [Hoerl and Kennard in 1970](https://doi.org/10.2307%2F1267351). It produces biased estimates of regression coefficients which are "shrunk" **towards zero** (but not equal to zero). Ridge regression finds the optimal regression coefficients `\(\hat\beta_{\text{ridge}}\)` by minimizing:

<div>
$$
\sum_{i=1}^n \left[ y_i - \left( \beta_0 + \sum_{j=1}^p \beta_j x_{ij} \right) \right]^2 + \lambda \sum_{j=1}^p \beta_j^2
$$
</div>

where `\(\lambda\)` is a constant weight term, and larger `\(\lambda\)` leads to more shrinkage. The value of `\(\lambda\)` is typically selected via cross validation.

From a Bayesian perspective, consider the following prior specification for the Bayesian linear regression model:

<div>
$$
p(\beta_0, \beta_1, \cdots, \beta_p, \sigma^2) = p(\beta_0) p(\beta_1, \cdots, \beta_p \mid \sigma^2) p(\sigma^2),
$$
</div>

where

`\begin{gathered}
\beta_0 \sim \mathcal{N}(0, \kappa^2) \text{ for some } \kappa > 0, \\
\beta_1, \cdots, \beta_p \mid \sigma^2 \overset{iid}{\sim} \mathcal{N} \left(0, \frac{\sigma^2}{\lambda} \right), \\
\sigma^2 \sim \text{IG}(0.01, 0.01)
\end{gathered}`

Hoerl and Kennard showed that the posterior modes for regression coefficients under this prior are equivalent to `\(\hat\beta_{\text{ridge}}\)`, the solution obtained via traditional ridge regression.

A benefit here is that we obtain **posteriors** for our biased regression parameters and could compute credible intervals, which is more challenging in the frequentist formulation. We do lose the pure interpretability of the coefficients themselves as we're introducing bias into those values.

Instead of fixing `\(\lambda\)` ahead of time (or use cross validation in the traditional setting), from the Bayesian perspective it might make more sense to place a prior distribution on `\(\lambda\)`, e.g.

`\begin{gathered}
\beta_1, \cdots, \beta_p \mid \sigma^2, \lambda \overset{iid}{\sim} \mathcal{N} \left(0, \frac{\sigma^2}{\lambda} \right), \\
\sigma^2 \sim \text{IG}(0.01, 0.01) \\
\lambda \sim \text{Gamma}(a, b)
\end{gathered}`

We may use uninformative priors such as `\(a=1\)` and `\(b=2\)` for `\(\lambda\)`.


#### LASSO

The Least Absolute Shrinkage and Selection Operator (LASSO) finds the optimal regression coefficients `\(\hat\beta_{\text{LASSO}}\)` by minimizing:

<div>
$$
\sum_{i=1}^n \left[ y_i - \left( \beta_0 + \sum_{j=1}^p \beta_j x_{ij} \right) \right]^2 + \lambda \sum_{j=1}^p |\beta_j|
$$
</div>

This method was popularized by [Tibshirani in 1996](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x). It simultaneously estimates regression coefficients while setting the ones that are "unimportant" **equal to zero**, i.e. it does variable selection. Similar to ridge regression, a larger `\(\lambda\)` leads to more regression coefficients being set to zero, and the value of `\(\lambda\)` is typically set via cross validation.

In Bayesian LASSO, consider the following priors:

<div>
$$
p(\beta_0, \beta_1, \cdots, \beta_p, \sigma^2) = p(\beta_0) p(\beta_1, \cdots, \beta_p \mid \sigma^2) p(\sigma^2),
$$
</div>

where

`\begin{gathered}
\beta_0 \sim \mathcal{N}(0, \kappa^2) \text{ for some } \kappa > 0 \\
\beta_1, \cdots, \beta_p \mid \sigma^2 \overset{iid}{\sim} \text{Laplace}\left( 0, \frac{\sqrt{\sigma^2}}{\lambda} \right) \\
\sigma^2 \sim \text{IG}(0.01, 0.01)
\end{gathered}`

[Park and Casella (2008)](https://doi.org/10.1198/016214508000000337) Showed that the posterior modes for regression coefficients under this prior are equivalent to `\(\hat\beta_{\text{LASSO}}\)`, the solution obtained via the traditional LASSO formulation.


### Some notes

From a Bayesian perspective, the posterior mode results are interesting but not necessarily useful, since one of the benefits of a Bayesian analysis is obtaining a *full posterior distribution*, rather than just a single point estimate.

In particular, this means that the Bayesian LASSO doesn't really perform variable selection, but it does still shrink regression parameters for improved fit. It also makes more sense to place a prior distribution on `\(\lambda\)` instead of fixing it, e.g. as suggested by Park and Casella:

`\begin{gathered}
\beta_1, \cdots, \beta_p \mid \sigma^2, \lambda \overset{iid}{\sim} \text{Laplace} \left( 0, \frac{\sqrt{\sigma^2}}{\lambda} \right) \\
\sigma^2 \sim \text{IG}(0.01, 0.01) \\
\lambda^2 \sim \text{Gamma}(a, b)
\end{gathered}`


### Diabetes example

Now that we have the theories explained, let's look at the posterior credible intervals for regression coefficients under the original prior specification, Bayesian ridge, and Bayesian LASSO of the diabetes example[^jags-code].



[^jags-code]: The R code for obtaining the posterior samples:

    
    ```r
    # Specify data, parameters, and initial values
    p <- ncol(X_train)
    dataList <- list(
      "n_train" = nrow(X_train),  # training sample size
      "p" = p,                    # number of variables
      "Y" = diab_train$y,         # response variable
      "X_train" = X_train,        # training data
      "n_test" = nrow(X_test),    # test sample size
      "X_test" = X_test           # test data
    )
    
    # JAGS settings
    adaptSteps <- 1000     # number of steps to "tune" the samplers
    burnInSteps <- 1000    # number of steps to "burn-in" the samplers
    nChains <- 3           # number of chains to run
    numSavedSteps <- 1000  # total number of steps in chains to save
    thinSteps <- 20        # number of steps to "thin" (1 = keep every step)
    nIter <- ceiling((numSavedSteps*thinSteps)/nChains)  # steps per chain
    
    run_jags_model <-
      function(model,
               init.vals,
               params,
               model.data = dataList,
               num.chains = nChains,
               adapt.steps = adaptSteps,
               burnin.steps = burnInSteps,
               num.iters = nIter,
               thin.steps = thinSteps,
               dic.mode = F) {
        # Define model
        m <- textConnection(model)
        jagsModel <- jags.model(
          m,
          data = model.data,
          inits = init.vals,
          n.chains = num.chains,
          n.adapt = adapt.steps
        )
        close(m)
        
        # Burn-in
        if (burnin.steps > 0) {
          update(jagsModel, n.iter = burnin.steps)
        }
        
        if (dic.mode) {
          # Generate DIC samples
          dic.samples(
            jagsModel, 
            variable.names = params, 
            n.iter = num.iters, 
            thin = thin.steps,
            type = "pD")
        } else {
          # Return posterior samples
          coda.samples(
            jagsModel,
            variable.names = params,
            n.iter = num.iters,
            thin = thin.steps
          )
        }
      }
    
    #########################
    # (1) Linear regression #
    #########################
    params1 <- c("beta_0", "beta", "sig2", "Y_pred")
    
    initValues1 <- list(
      "beta_0" = 0,
      "beta" = rep(0, p),
      "tau2" = 1
    )
    
    # Run model in JAGS
    trad_diab <- "
    model {
      # Sampling model
      for (i in 1:n_train) {
      	Y[i] ~ dnorm(mu[i], tau2)
        mu[i] = beta_0 + inprod(beta, X_train[i, ])
      }
    
      # Prior model
      beta_0 ~ dnorm(0, 1e-10)
      for (j in 1:p) {
        beta[j] ~ dnorm(0, 1)
      }
      tau2 ~ dgamma(0.01, 0.01)
    
      # Need to have model calculate sig2 = 1/tau2
      sig2 = 1/tau2
    
      # Prediction on test data
      for (i in 1:n_test) {
    	  mu_pred[i] = beta_0 + inprod(beta, X_test[i, ])
    	  Y_pred[i] ~ dnorm(mu_pred[i], tau2)
      }
    }
    "
    
    codaSamples1 <- run_jags_model(trad_diab, initValues1, params1)
    
    ########################
    ## (2) Bayesian ridge ##
    ########################
    params2 <- c("beta_0", "beta", "sig2", "lambda", "Y_pred")
    
    initValues2 <- list(
      "beta_0" = 0,
      "beta" = rep(0, p),
      "tau2" = 1,
      "lambda" = 0.1
    )
    
    # Run model in JAGS
    ridge_diab <- "
    model {
      # Sampling model
      for (i in 1:n_train) {
      	Y[i] ~ dnorm(mu[i], tau2)
        mu[i] = beta_0 + inprod(beta, X_train[i, ])
      }
    
      # Prior model
      beta_0 ~ dnorm(0, 1e-10)
    
      for (j in 1:p) {
        beta[j] ~ dnorm(0, lambda*tau2)
      }
    
      lambda ~ dgamma(1, 2)
      tau2 ~ dgamma(0.01, 0.01)
    
      # Need to have model calculate sig2 = 1/tau2
      sig2 = 1/tau2
    
      # Prediction on test data
      for (i in 1:n_test) {
    	mu_pred[i] = beta_0 + inprod(beta, X_test[i, ])
    	Y_pred[i] ~ dnorm(mu_pred[i], tau2)
      }  
    }
    "
    
    codaSamples2 <- run_jags_model(ridge_diab, initValues2, params2)
    
    ########################
    ## (3) Bayesian LASSO ##
    ########################
    params3 <- c("beta_0", "beta", "sig2", "lambda", "Y_pred")
    
    initValues3 <- list(
      "beta_0" = 0,
      "beta" = rep(0, p),
      "tau2" = 1,
      "lambda2" = 0.1
    )
    
    # Run model in JAGS
    lasso_diab <- "
    model {
      # Sampling model
      for (i in 1:n_train) {
      	Y[i] ~ dnorm(mu[i], tau2)
        mu[i] = beta_0 + inprod(beta, X_train[i, ])
      }
    
      # Prior model
      beta_0 ~ dnorm(0, 1e-10)
    
      for(j in 1:p){
        beta[j] ~ ddexp(0, sqrt(lambda2*tau2))
      }
    
      lambda2 ~ dgamma(1,2)
      tau2 ~ dgamma(0.01,0.01)
    
      # Need to have model calculate sig2 and lambda
      sig2 = 1/tau2
      lambda = sqrt(lambda2)
    
      # Prediction on test data
      for (i in 1:n_test) {
    	mu_pred[i] = beta_0 + inprod(beta, X_test[i, ])
    	Y_pred[i] ~ dnorm(mu_pred[i], tau2)
      }
    }
    "
    codaSamples3 <- run_jags_model(lasso_diab, initValues3, params3)
    ```


```r
par(mar = c(2, 1, 2, 1), mfrow = c(1, 3))
MCMCplot(codaSamples1, params = "beta", main = "(a)",
         ref = NULL, labels = NULL, xlim = c(-2, 2))
MCMCplot(codaSamples2, params = "beta", main = "(b)",
         ref = NULL, labels = NULL, xlim = c(-2, 2))
MCMCplot(codaSamples3, params = "beta", main = "(c)",
         ref = NULL, labels = NULL, xlim = c(-2, 2))
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/caterpillar-plots-1.svg" alt="Posterior CIs for beta under (a) original prior, (b) Bayesian ridge, and (c) Bayesian LASSO." width="1440" />
<p class="caption">Figure 2: Posterior CIs for beta under (a) original prior, (b) Bayesian ridge, and (c) Bayesian LASSO.</p>
</div>

In Figure <a href="#fig:caterpillar-plots">2</a>(a) we have the posterior CIs under the traditional Bayesian linear regression setting. Note that the first ten intervals correspond to the main effects, and the rest are interaction and quadratic terms. Typically with those regression models, we hope many of the higher-order terms to be close to zero, because we probably don't expect them to explain `\(y\)` much.

We can see that there's a lot of variability in the coefficient values. One of the reasons for this is we have a large number of predictor variables relative to the sample size. Under the Bayesian ridge setting as shown in Figure <a href="#fig:caterpillar-plots">2</a>(b), most of the coefficients are pushed towards zero. The posterior CIs are also narrower for a lot of the coefficients. In Bayesian LASSO we observe similar effects.

The coefficients are more concentrated around zero with less variability due to shrinkage. We're essentially specifying a prior saying the coefficients are centered at zero and have values similar to each other, meaning that their values are all going to be similar to zero. This is a stronger prior than the previous model.


#### Comparison of posteriors

So what effect does this have on individual regression coefficients? When the number of predictors is large, the estimation of the coefficients can be unstable. This can be seen by looking at trace plots and posterior densities of the regression coefficients. Here we take `\(\beta_{57}\)` as an example.


```r
ggMCMCtrace <- function(mcmc.chains, param, num.chains) {
  # ggplot version of MCMCvis::MCMCtrace()
  library(patchwork)
  
  num.iters <- nrow(mcmc.chains) %/% num.chains
  dat <- mcmc.chains %>%
    as_tibble() %>%
    select(!!param) %>%
    mutate(
      Chain = factor(rep(seq(num.chains), each = num.iters)),
      Iteration = rep(seq(num.iters), times = num.chains)
    )
  
  p1 <- ggline(dat, x = "Iteration", y = param, color = "Chain",
               palette = "nejm", plot_type = "l", numeric.x.axis = T)
  p2 <- ggdensity(dat, x = param,
                  xlab = "Parameter estimate", ylab = "Density")
  
  p1 | p2
}

mcmcChain1 <- as.matrix(codaSamples1)
ggMCMCtrace(mcmcChain1, param = "beta[57]", num.chains = nChains)
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/beta57-trace-density-1.svg" alt="Trace plot and posterior density of a regression coefficient under the original prior." width="960" />
<p class="caption">Figure 3: Trace plot and posterior density of a regression coefficient under the original prior.</p>
</div>

The trace plot doesn't look great -- the instability in the chains comes from increased variation in estimating regression coefficients when `\(p\)` is large. The trace plots of the Bayesian ridge and LASSO show much more of the convergent behavior because we are restricting the values of `\(\beta\)` by introducing bias.


```r
mcmcChain2 <- as.matrix(codaSamples2)
mcmcChain3 <- as.matrix(codaSamples3)

ggMCMCtrace(mcmcChain2, param = "beta[57]", num.chains = nChains) /
  ggMCMCtrace(mcmcChain3, param = "beta[57]", num.chains = nChains)
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/beta57-ridge-lasso-1.svg" alt="Trace plot and posterior density of a regression coefficient under the Bayesian ridge (top) and Bayesian LASSO (bottom)." width="960" />
<p class="caption">Figure 4: Trace plot and posterior density of a regression coefficient under the Bayesian ridge (top) and Bayesian LASSO (bottom).</p>
</div>

For this specific coefficient, we have from the density plots `\(\hat\beta_{\text{OLS}, 57} = 0.4\)` and `\(\hat\beta_{\text{ridge}, 57} = \hat\beta_{\text{LASSO}, 57} = 0.15\)`. The ridge and LASSO gave different posteriors for the regression coefficient, but both are pushed towards zero as we've seen in the CI plots.


## Model comparison

One traditional way to compare these three models is to quantify how well they predict disease progression response on a test sample. We can predict the response value `\(\hat{y}_i\)` under each model, and compare to the true value `\(y_i\)`.


### Posterior-based predictions

For the `\(i\)`-th patient in the test sample, we have its vector of covariantes `\(x_i\)`. The `posterior predictive distribution` for its response is:

<div>
$$
\boldsymbol{Y}_{\text{pred, }i} \mid \boldsymbol{\beta}, \sigma^2 \sim \mathcal{N}(\boldsymbol{x}_i^T \boldsymbol{\beta}, \sigma^2)
$$
</div>

We sampled this in JAGS and stored the values in the `Y_pred` column. We get a whole distribution for the predictions for each patient, and we may consider taking `\(\hat{y}_i\)` to be the `posterior predictive mean`, i.e.

<div>
$$
\hat{y}_i = E\left[ \boldsymbol{Y}_{\text{pred, } i} \mid \boldsymbol{\beta}, \sigma^2 \right]
$$
</div>

We can then quantify how well each of the three models predict our test data by computing the `squared-error loss`:

<div>
$$
\text{prediction error} = \frac{1}{n_{\text{test}}} \sum_{i=1}^{n_\text{test}} (\hat{y}_i - y_i)^2
$$
</div>

The prediction errors for the three models are summarized in the table below. Both of the penalized regression methods demonstrate improved predictive performance on the test samples[^ridge-vs-lasso].


|Prior      | Prediction error|
|:----------|----------------:|
|Original   |           0.6057|
|Ridge-type |           0.5467|
|LASSO-type |           0.5320|

[^ridge-vs-lasso]: In this case the LASSO seems to perform better than the ridge, but often times the ridge performs slightly better. The LASSO is usually more interpretable as it sets some coefficients to zero.


### Model selection

We compared the three models above by how good their predictions were. Given several candidate models, how can we select the "best" model? In general, we tend to favor models which are:

- Scientifically reasonable, perhaps backed by existing theory.
- Performs well at predicting test data. Sometimes this is not the model that makes most scientific sense!
- Interpretable.
- Parsimonious, i.e. we generally prefer simpler models.

Different applications may place different weights on these criterion, and as a result, there are many ways to evaluate and select models. Here we introduce a few popular Bayesian techniques for `model selection`. Justifying the use of a model may involve using one or more of those techniques in conjunction with other considerations that are application-dependent.


#### Bayes factor

Suppose we have two competing models denoted `\(M_1\)` and `\(M_2\)`, with respective parameter vectors `\(\theta_1\)` and `\(\theta_2\)`. The `Bayes factor` of `\(M_1\)` to `\(M_2\)` is defined as:

<div>
$$
\text{BF} = \frac{\text{posterior odds of model 1 to model 2}}{\text{prior odds of model 1 to model 2}}
$$
</div>

It can be shown that the above formula reduces to the following:

<div>
$$
\text{BF} = \frac{p(y \mid M_1)}{p(y \mid M_2)}
$$
</div>

where

<div>
$$
p(y \mid M_i) = \int p(y \mid \theta_i, M_i) p(\theta_i \mid M_i) d\theta_i
$$
</div>

is the `marginal likelihood` of our data under model `\(M_i\)`[^model-evidence]. In the integral we have the sampling model for model `\(M_i\)` and the prior for `\(\theta_i\)` under `\(M_i\)`.

[^model-evidence]: In Bayes' Theorem:

    <div>
    $$
    p(\theta \mid y) = \frac{p(y \mid \theta) p(\theta)}{p(y)}
    $$
    </div>
    
    We often ignore the denominator part, which is exactly this marginal likelihood value. This is sometimes called the model `evidence`.

The idea is if the Bayes factor is greater than one, then we would prefer model 1 to model 2 as the marginal likelihood of data under model 1 is higher. A rule of thumb is `\(\text{BF} > 100\)` is strong evidence in favor of model `\(M_1\)`.

A few caveats with Bayes factors:

- Bayes factor only compares two models.
- The marginal likelihoods computed in the Bayes factor are generally very difficult to compute or approximate unless priors are conjugate. 
- Bayes factors cannot be used for models with improper priors as the marginal likelihood can't be computed.

Bayes factors are also problematic if models are nested with the dimension of `\(\theta_1 < \theta_2\)`, e.g.

`\begin{gathered}
M_1: \mu_i = \beta_0 + \beta_1 x_{1i} \\
M_2: \mu_i = \beta_0 + \beta_1 x_{1i} + \beta_2 x_{2i}
\end{gathered}`

`\(M_1\)` is nested within `\(M_2\)`, since we could set `\(\beta_2 = 0\)` in `\(M_2\)` to get `\(M_1\)`. If we have a large number of predictors, we often want to see whether including a predictor would increase model performance.


#### Deviance information criterion

Another thing we might do is to compute model-based criterion measures (similar to AIC and BIC), which systematically assign scores to different models. The score should penalize models with a large number of parameters. This is generally beneficial when models get complex and methods like stepwise regression becomes difficult.

Suppose we are interested in estimating parameters `\(\theta\)` under a model with sampling density `\(p(y \mid \theta)\)`. Let `\(\hat\theta\)` be the `posterior mean`, perhaps estimated using MCMC. The `deviance information criterion` (DIC) for this model is defined as:

<div>
$$
\text{DIC} = -2\log \left( p(y \mid \hat\theta) \right) + 2p_D
$$
</div>

where the first term is called the `deviance`, and the second term is an estimate of the **effective number of model parameters**, which is used to penalize against overly complex models. Models with **smaller** values of DIC are preferred.

The DIC is preferred from a Bayesian perspective compared to other criterion like AIC or BIC because:

- DIC depends on the posterior mean, which incorporates prior information.
- DIC can be easily computed from samples generated via MCMC.
- AIC and BIC are more "frequentist-driven", as they both rely on maximum likelihood estimates instead of the posterior mean.

In JAGS, the function `dic.samples` generates samples from the posterior while simultaneously computing the DIC for the model[^diabetes-dic].

[^diabetes-dic]: The DICs of the three models in the diabetes example are:

    
    ```
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 342
    ##    Unobserved stochastic nodes: 166
    ##    Total graph size: 30131
    ## 
    ## Initializing model
    ```
    
    ```
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 342
    ##    Unobserved stochastic nodes: 167
    ##    Total graph size: 30134
    ## 
    ## Initializing model
    ```
    
    ```
    ## Compiling model graph
    ##    Resolving undeclared variables
    ##    Allocating nodes
    ## Graph information:
    ##    Observed stochastic nodes: 342
    ##    Unobserved stochastic nodes: 167
    ##    Total graph size: 30136
    ## 
    ## Initializing model
    ```
    
    
    
    | Original|  Ridge|  LASSO|
    |--------:|------:|------:|
    |   779.49| 763.42| 761.33|

#### Other methods

**Stochastic search methods** are often used for large model spaces, i.e. many candidate models arising from the large number of predictors. **Predictive checking** looks at posterior predictive distributions and compare those to true values. Cross validation is an example of this.
