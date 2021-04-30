---
title: "Bayesian Generalized Linear Models"
date: 2021-04-19T10:20:00-05:00
summary: 'This lecture discusses a simple logistic regression model for predicting a binary variable. GLMs are necessary when the response variable cannot be modeled appropriately by a normal distribution, and use a link function to connect parameters of the response distribution to the covariates.' # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-generalized-linear-models"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 120 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Sparrow with a twig." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "njOwWb0v6Vs" # Unsplash ID of the picture
---



This lecture briefly introduces how to fit a generalized linear model from a Bayesian perspective. Prior specification proceeds similar to Bayesian linear regression, and JAGS can be used fairly simply to sample from the posterior.

Although we will be discussing a fairly simple example, the GLM framework could be combined with some aforementioned topics such as using random intercepts to aggregate data at a group level.

## Sparrow example

This motivating example comes from Hoff's book, page 244. Younger male sparrows may or may not nest during a mating season, based on their physical characteristics. Researchers have recorded the nesting success of 43 young male sparrows of the same age, as well as their wingspan in centimeters.


```r
library(tidyverse)
library(ggpubr)
library(rjags)
```


```r
sparrows <- read_table(
  "http://www2.stat.duke.edu/courses/Fall09/sta290/datasets/Hoffdata/msparrownest.dat",
  col_names = c("success", "ws")
)
```

The goal is to model the relationship between nesting success and wingspan. Our response data is binary:

<div>
$$
    Y_i = \begin{cases}
        1, & \text{if sparrow }i \text{ successfully nested,} \\
        0, & \text{otherwise}
    \end{cases}
$$
</div>

It seems we tend to observe greater wingspans among sparrows which successfully nested.


```r
ggboxplot(sparrows, x = "success", y = "ws",
          xlab = "Nesting success", ylab = "Wingspan (cm)")
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/boxplot-sparrow-1.svg" alt="Boxplot of wingspan by indicator of nesting success." width="480" />
<p class="caption">Figure 1: Boxplot of wingspan by indicator of nesting success.</p>
</div>

The response is no longer a continuous variable, so the normal distribution doesn't seem appropriate here.


## Logistic regression model

We observe the *binary* response variable `\(Y_i \in \{0, 1\}\)`, indicating whether sparrow `\(i\)` has successfully nested. For the **sampling model**, we can think of this as a Binomial experiment:

<div>
$$
    Y_i \mid \theta_i \sim \text{Ber}(\theta_i) = \text{Bin}(n=1, \theta_i) \text{ for } i = 1, \cdots, n
$$
</div>

where `\(\theta_i\)` is the probability of nesting success for sparrow `\(i\)` and is constrained to be between 0 and 1. We are going to assume conditional independence too, i.e. conditional independence of the responses based on the probabilities of success. The probability density function is thus:

<div>
$$
    p(y \mid \theta) = \prod_{i=1}^n p(y_i \mid \theta_i)
$$
</div>

The question is how do we relate this nesting success to wingspan, the covariate in this model? One idea is to treat this as a linear regression model, e.g. `\(Y_i \mid \theta_i \sim \text{Ber}(\theta_i)\)`, where

`\begin{equation}
  \theta_i = \beta_1 + \beta_2 x_i
  \tag{1}
\end{equation}`


The problem is `\(\theta_i\)` is a probability, so it must be **between 0 and 1**. The way the regression function is specified in Eq.(1) places no constraints for `\(\beta_1\)`, `\(\beta_2\)`, or a positive `\(x_i\)`. There's no guarantee that for a sparrow with wingspan `\(x_i\)`, we get a valid nesting probability `\(\theta_i\)` between 0 and 1.

From a Bayesian perspective, we may construct a prior distribution on `\(\beta_1\)` and `\(\beta_2\)` so that `\(\theta_i\)` *has* to be between 0 and 1, but it's not really a rational way to address this issue.

We may consider resolving this issue by mapping `\(\theta_i\)` to the real number line, i.e. from `\(-\infty\)` to `\(\infty\)`, through the following transformation:

<div>
$$
  g(\theta_i) = \text{logit}(\theta_i) = \log \left( \frac{\theta_i}{1 - \theta_i} \right)
$$
</div>

Here `\(g\)` is the `logit` or `log-odds` function. Note that:

`\begin{gathered}
  \theta_i \approx 0 \implies g(\theta_i) \rightarrow -\infty \\
  \theta_i \approx 1 \implies g(\theta_i) \rightarrow \infty
\end{gathered}`

`\(g\)` maps probabilities `\(\theta_i\)` to the real number line, which is exactly what we wanted! The issue could then be fixed by:

<div>
$$
  g(\theta_i) = \log \left( \frac{\theta_i}{1 - \theta_i} \right) = \beta_1 + \beta_2 x_i
$$
</div>

By doing this, we no longer need a constraint on `\(\beta_1\)`, `\(\beta_2\)`, or `\(x_i\)`. As for **prior models**, we don't need priors on `\(\theta_i\)`'s anymore, but we still need priors on `\(\beta_1\)` and `\(\beta_2\)`. Assuming we don't have much information:

<div>
$$
  \beta_1, \beta_2 \overset{iid}{\sim} \mathcal{N}(0, 10^{10})
$$
</div>


## Generalized linear models

The previous section is an example of what we call `generalized linear models`. GLMs are used to relate a linear function of the predictor variables to a response variable `\(Y\)` which is not normally distributed. The basic idea is that we specify a sampling model for the response variable `\(Y_i\)` with parameters `\(\eta_i\)`. We use a `link function` `\(g\)` to connect the parameters to the `\(p\)` predictors, i.e.,

<div>
$$
  g(\eta_i) = \beta_1 + \beta_2 x_{i, 2} + \cdots + \beta_{p+1}x_{i, p+1} = \beta_1 + \sum_{j=2}^{p+1} \beta_j x_{ij}
$$
</div>

The motivation behind this is that the parameter space of `\(\eta_i\)` might not be all real numbers, but it's possible for the RHS to take any real number. So `\(g\)` (and its inverse[^invertible-g]) maps this to the appropriate parameter space.

[^invertible-g]: One requirement of link functions is they should be invertible.

There are a lot of examples for GLMs, and two of the most common are `logistic regression` and `Poisson regression`. Logistic regression is used for binary (0-1) response data, and it uses a logit link. Poisson regression is used for count response data, and a log link is used:

`\begin{gathered}
  Y_i \mid \theta_i \sim \text{Pois}(\theta_i) \\
  g(\theta_i) = \log(\theta_i) = \beta_1 + \sum_{j=2}^{p+1} \beta_j x_{ij}
\end{gathered}`

The only constraint here is `\(\theta_i > 0\)`, as it represents the mean number of events occurring.


## Back to sparrows

The full Bayesian logistic regression model for the sparrow data is given below. For `\(i = 1, \cdots, n\)`, the **sampling model** is given by:

`\begin{gathered}
  Y_i \mid \theta_i \sim \text{Ber}(\theta_i) \\
  \log\left(\frac{\theta_i}{1-\theta_i}  \right) = \beta_1 + \beta_2 x_i
\end{gathered}`

The **prior model** is:

<div>
$$
  \beta_1, \beta_2 \overset{iid}{\sim} \mathcal{N}(0, 10^{10})
$$
</div>

We can draw approximate samples from the posterior using MCMC. Fitting GLMs is fairly simple in JAGS, as shown in the code below.


```r
n <- nrow(sparrows)  # sample size

## Logistic regression
# Specify data, parameters, and initial values
dataList <- list(
  "n" = n,
  "success" = sparrows$success,
  "ws" = sparrows$ws
)

parameters <- c("theta_pred", "beta_1", "beta_2")

initValues <- list(
  "beta_1" = 0,
  "beta_2" = 1
)

# JAGS settings
adaptSteps <- 10000     # number of steps to "tune" the samplers
burnInSteps <- 10000    # number of steps to "burn-in" the samplers
nChains <- 3            # number of chains to run
numSavedSteps <- 10000  # total number of steps in chains to save
thinSteps <- 100        # number of steps to "thin" (1 = keep every step)
nIter <- ceiling((numSavedSteps*thinSteps)/nChains) 	# steps per chain

# Run model in JAGS
m <- textConnection("
model {
  # Sampling model
  for (i in 1:n) {
  	success[i] ~ dbern(theta[i])
    logit(theta[i]) = beta_1 + beta_2*ws[i]
  }

  # Prior model
  beta_1 ~ dnorm(0, 1e-10)
  beta_2 ~ dnorm(0, 1e-10)

  # Prediction - suppose a sparrow has wingspan of 12.6
  ws_pred = 13.7
  logit(theta_pred) = beta_1 + beta_2*ws_pred
}
")
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
mcmcChain <- as.matrix(codaSamples)
```

Numerical summaries of the posterior regression coefficients are given below.


|Parameter | Posterior mean| Posterior sd|95% posterior CI |
|:---------|--------------:|------------:|:----------------|
|$\beta_1$ |         -11.01|         4.91|(-21.37, -2.03)  |
|$\beta_2$ |           0.87|         0.38|(0.17, 1.68)     |

Note that the 95% credible interval for `\(\beta_2\)` doesn't include 0, which is strong evidence suggesting that wingspan is related to nesting success.

Interpretation of the coefficients is still possible, just somewhat more complicated than linear regression. `\(E(\beta_1 \mid y) = -11.01\)` means that for wingspan of `\(x_i = 0\)`, the log-odds of nesting success is -11.01:

<div>
$$
  g(\theta_i) = \log\left( \frac{\theta_i}{1 - \theta_i} \right) = -11.01 \implies \theta_i = \frac{\exp(-11.01)}{1 + \exp(-11.01)} \approx 0
$$
</div>


`\(E(\beta_2 \mid y) = 0.87\)` means for an increase in wingspan of 1cm, we expect the log-odds of nesting success to increase by 0.87. The odds is expected to increase by a factor of `\(\exp(0.87) = 2.39\)`.


### Predictions

Suppose we have a new sparrow with a wingspan of 13.7cm. Do we expect this sparrow to be able to successfully nest? Given `\(x^* = 13.7\)`, one way to answer this question is by looking at its posterior probability of nesting success, based on the observed data.

<div>
$$
  \log \left( \frac{\theta_{\text{pred}}}{1 - \theta_{\text{pred}}} \right) = \beta_1 + \beta_2 x^*
$$
</div>

We may get `\(\theta_{\text{pred}}\)` from the posteriors of `\(\beta_1\)` and `\(\beta_2\)`. The posterior mean is 0.71 and the 95% CI is (0.52, 0.87), which suggest it's more likely than not this particular sparrow is going to nest.


### Notes

We can build GLMs with multiple predictors, and use similar ideas from Bayesian linear regression to think about model selection, i.e. DIC. We can account for **random effects** by combining GLMs with ideas from hierarchical modeling. Keeping the sampling model unchanged, the new regression function is:

<div>
$$
  \log \left( \frac{\theta_i}{1 - \theta_i} \right) = \alpha_{s_i} + \beta_1 + \beta_2 x_i
$$
</div>


where `\(s_i\)` represents the nesting location of sparrow `\(i\)`, and `\(\alpha_{s_i}\)` are the random intercepts. The "location-level model" is:

<div>
$$
  \alpha_1, \cdots, \alpha_J \mid \sigma_\alpha^2 \overset{iid}{\sim} \mathcal{N}(0, \sigma_\alpha^2)
$$
</div>

Priors for `\(\beta_1\)` and `\(\beta_2\)` are also the same as before.

There are entire courses devoted to Bayesian GLMs, and this lecture is only meant to be a simple introduction of this topic.
