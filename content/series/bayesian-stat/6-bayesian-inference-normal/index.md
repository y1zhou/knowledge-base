---
title: "Bayesian Inference for the Normal Model"
date: 2021-02-22T10:20:00-05:00
summary: "The normal distribution has two parameters, but we focus on the one-parameter setting in this lecture. We also introduce the posterior predictive check as a way to assess model fit, and briefly discuss the issue with improper prior distributions." # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-normal-distribution"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 60 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Dominoes (more importantly, numbers)." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "ePHz9WOME0c" # Unsplash ID of the picture
---

This lecture discusses Bayesian inference of the normal model, most commonly used for continuous data. The normal distribution has two parameters (corresponding to mean and variance), and so while we will ultimately discuss posterior inference for multiple parameters, we focus on posterior estimation of a single parameter in this lecture -- in particular, estimating the mean given a fixed variance.

The example we're going to study is about crab shells. Two scientists[^campbell-mahon] were interesting in studying mophological variation of _Leptograpsus variegatus_, also known as the purple rock crab. They collected $n=200$ crabs and measured their shape through 5 variables. Our focus will be on **carapace length** (in mm), the length of the crab's protective upper shell. The goal is to estimate the typocal carapace length, as well as the variability in lengths.

[^campbell-mahon]: Campbell, N.A. and Mahon, R.J., 1974. A multivariate study of variation in two species of rock crab of the genus Leptograpsus. _Australian Journal of Zoology_, 22(3), pp.417-425.

## Exploratory data analysis

The dataset can be obtained from the `MASS` R package. The `crabs` data frame has 200 rows and 8 columns, describing 5 morphological measurements on 50 crabs each of two color forms and both sexes.

From the histogram of the carapace lengths[^crab-hist], we can see the values range from 10 to 50 millimeters and is centered around 30mm. The mean is 32.11mm and the variance is 50.68mm$^2$. The distribution looks overall symmetric with no outliers.

{{< figure src="crab_length_hist.png" caption="Carapace length of *Leptograpsus variegatus*." numbered="true" >}}

[^crab-hist]: R code for plotting:

    ```r
    library(ggpubr)

    data(crabs, package = "MASS")
    mean(crabs$CL)
    var(crabs$CL)

    gghistogram(crabs, x = "CL", add_density = T,
                xlab = "Carapace length (mm)", ylab = "Frequency")
    ```

## Sampling model

Let $Y$ be a random variable representing the carapace length (for one of the crabs). The `normal model` is often used for continuous quantities which are _distributed symmetrically_ with _few outliers_:

<div>
$$
    Y \mid \theta, \sigma^2 \sim N(\theta, \sigma^2)
$$
</div>

The parameters of interest are $\theta$, the true average carapace length (in mm); and $\sigma^2$, the true variance in carapace length (in mm$^2$). The probability density function is:

<div>
$$
    p(y \mid \theta, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp \left\{ -\frac{(y-\theta)^2}{2\sigma^2} \right\}, \quad -\infty < y < \infty
$$
</div>

Like the Poisson model, our data actually consists of $n=200$ observations, each representing the carapace length of a different crab. We will assume that each observation is **conditionally independent**:

<div>
$$
    Y_1, \cdots, Y_n \mid \theta, \sigma^2 \overset{iid}{\sim} N(\theta, \sigma^2)
$$
</div>

This means if `$\boldsymbol{Y} = (Y_1, \cdots, Y_n)$`, then

<div>
$$
    \begin{aligned}
        p(Y \mid \theta, \sigma^2) &= \prod_{i=1}^n p(Y_i \mid \theta, \sigma^2) \\
        &= \prod_{i=1}^n \frac{1}{\sqrt{2\pi\sigma^2}} \exp \left\{ -\frac{(Y_i - \theta)^2}{2\sigma^2} \right\} \\
        &= (2\pi\sigma^2)^{-\frac{n}{2}} \exp \left\{ -\frac{1}{2\sigma^2} \sum_{i=1}^n (Y_i - \theta)^2 \right\}
    \end{aligned}
$$
</div>

### Estimation of multiple parameters

This is the first model where we want to estimate multiple parameters. There are two ways to approach inference:

1. Assume one parameter is _known_, and compute the posterior for the other parameter. This is an overly simplistic view, but this is our focus in this lecture.

    When $\theta$ is unknown and $\sigma^2$ is known, we fix $\sigma^2$:

    <div>
    $$
        p(\theta \mid y, \sigma^2) \propto p(y \mid \theta, \sigma^2) p(\theta)
    $$
    </div>

    When $\theta$ is known and $\sigma^2$ is unknown:

    <div>
    $$
        p(\sigma^2 \mid y, \theta) \propto p(y \mid \theta, \sigma^2) p(\sigma^2)
    $$
    </div>

2. Assume both parameters are _unknown_, and compute a joint posterior for $(\theta, \sigma^2)$:

    <div>
    $$
        p(\theta, \sigma^2 \mid y) \propto p(y \mid \theta, \sigma^2) p(\theta, \sigma^2)
    $$
    </div>

    We will discuss this two-dimensional prior case in later lectures.

## Prior model for $\theta$

For the rest of this lecture, we assume that the scientist is extremely confident that the true variance $\sigma^2 = 45 mm^2$, i.e. we will treat this value as _fixed and known_. In this case, we only need to specify a prior model for $\theta$, the average carapace length. One suitable choice is:

<div>
$$
    \theta \sim N(\mu_0, \tau_0^2)
$$
</div>

where `$\mu_0$` is chosen as the "center" of the prior beliefs, and `$\tau_0^2$` reflects uncertainty in prior beliefs. As we increase `$\tau_0^2$`, the distribution becomes more spread out. The probability density function for the prior distribution is:

<div>
$$
    p(\theta) = \frac{1}{\sqrt{2\pi\tau_0^2}} \exp \left\{ -\frac{(\theta - \mu_0)^2}{2\tau_0^2} \right\}
$$
</div>

By Bayes' Theorem, the posterior can be computed:

<div>
$$
    \begin{aligned}
        p(\theta \mid y, \sigma^2) &= p(y \mid \theta, \sigma^2) p(\theta) \\
        &= p(\theta) p(y_1 \mid \theta, \sigma^2) \cdots p(y_n \mid \theta, \sigma^2) \\
        &= \frac{1}{\sqrt{2\pi\tau_0^2}} \exp \left\{ -\frac{(\theta - \mu_0)^2}{2\tau_0^2} \right\} \times (2\pi\sigma^2)^{-\frac{n}{2}} \exp \left\{ -\frac{1}{2\sigma^2} \sum_{i=1}^n (y_i - \theta)^2 \right\} \\
        &\propto \exp \left\{ -\frac{(\theta - \mu_0)^2}{2\tau_0^2} \right\} \times \exp \left\{ -\frac{1}{2\sigma^2} n(\bar{y} - \theta)^2 - \frac{1}{2\sigma^2} \sum_{i=1}^n (y_i - \bar{y})^2 \right\} \\
        &\propto \exp \left\{ -\frac{(\theta - \mu_0)^2}{2\tau_0^2} \right\} \times \exp \left\{ -\frac{n(\bar{y} - \theta)^2}{2\sigma^2} \right\} \\
        &\propto \exp \left\{ -\frac{\theta^2 - 2\theta\mu_0 + \mu_0^2}{2\tau_0^2} - \frac{n\bar{y}^2 - 2n\bar{y}\theta + n\theta^2}{2\sigma^2} \right\} \\
        &\propto \exp \left\{ -\frac{\theta^2}{2}\left( \frac{1}{\tau_0^2} + \frac{n}{\sigma^2} \right) + \theta \left( \frac{\mu_0}{\tau_0^2} + \frac{n\bar{y}}{\sigma^2} \right) - \frac{1}{2}\left( \frac{\mu_0^2}{\tau_0^2} + \frac{n\bar{y}^2}{\sigma^2} \right) \right\}
    \end{aligned}
$$
</div>

To make this into a kernel of a normal distribution, we want it to be in the form:

<div>
$$
    \exp\left\{ -\frac{1}{2\tau_n^2} (\theta - \mu_n)^2 \right\}
    = \exp\left\{ -\frac{\theta^2}{2} \cdot \frac{1}{\tau_n^2} + \theta \cdot \frac{\mu_n}{\tau_n^2} -\frac{1}{2} \cdot \frac{\mu_n^2}{\tau_n^2} \right\}
$$
</div>

We can find that the posterior follows $N(\mu_n, \tau_n^2)$, where

<div>
$$
    \mu_n = \cfrac{\cfrac{\mu_0}{\tau_0^2} + \cfrac{n\bar{y}}{\sigma^2}}{\cfrac{1}{\tau_0^2} + \cfrac{n}{\sigma^2}}, \quad \tau_n^2 = \cfrac{1}{\cfrac{1}{\tau_0^2} + \cfrac{n}{\sigma^2}}
$$
</div>

and $\mu_n$ is found from its ratio with $\tau_n^2$ in the second term.

### Prior hyperparameter selection for $\theta$

Suppose we are very confident based on previous studies that the average carapace length is somewhere between 12 and 44 mm. How can we select the prior hyperparameters $\mu_0$ and $\tau_0^2$?

With this prior belief, it makes sense to set the center at the middle of 12 and 44 mm, and that midpoint is 28mm. Now the question is how much uncertainty should be associated with the prior here. If we are 95% confident, then going at 2 standard deviations from $\mu_0$ gives us 95% of the probability:

<div>
$$
    \tau_0 = s.d. = \frac{44 - 28}{2} = 8
$$
</div>

Thus the prior is set to be $\theta \sim N(28, 8^2)$. Note that prior information is very subjective, and it's not always easy to interpret the hyperparameters. If we are, say, 99% confident, then we may go out 3 sd's instead of 2. Ultimately because we have 200 data points, this shouldn't affect our posterior too much as it's dominated by the observed data. Below is a plot of $S=1000$ Monte Carlo samples from our prior distribution generated via:

<div>
$$
    \theta^{(s)} \sim N(\mu_0 = 28, \tau_0^2 = 64), \quad s = 1, 2, \cdots, S
$$
</div>

{{< figure src="MC_prior.png" caption="Histogram of 1000 Monte Carlo samples generated from the prior." numbered="true" >}}

## Posterior and posterior summaries

Using the formula for $\mu_n$ and $\tau_n^2$, the posterior mean and variance are

<div>
$$
    \begin{gathered}
        E(\theta \mid y, \sigma^2) = 32.1mm \\
        Var(\theta \mid y, \sigma^2) = 0.22 mm^2
    \end{gathered}
$$
</div>

The 95% credible interval for $\theta$ is (31.18mm, 33.06mm). Although we already know the true posterior distribution, we can use MC samples to find the mean and variance.

### Posterior mean and variance decomposition

Just like in previous models, we can write the posterior mean for $\theta$ as a weighted average of the _prior mean_ and _sample mean_. We can also find a similar decomposition for a quantity related to the posterior variance for $\theta$. The formulas are messy, and they become much easier to write in terms of the `sampling and prior precisions`:

<div>
$$
    \begin{gathered}
        \text{sampling precision}: \tilde{\sigma}^2 = \frac{1}{\sigma^2} \\
        \text{prior precision}: \tilde{\tau}_0^2 = \frac{1}{\tau_0^2}
    \end{gathered}
$$
</div>

Precision is just defined as the inverse of the variance. Whereas variance is interpreted as "degree of variability", precision can be thought of as "amount of information". With these quantities, the `posterior mean decomposition` is

<div>
$$
    \mu_n = \frac{\tilde{\tau}_0^2}{\tilde{\tau}_0^2 + n \tilde{\sigma}^2} \mu_0 + \frac{n\tilde{\sigma}^2}{\tilde{\tau}_0^2 + n\tilde{\sigma}^2} \bar{y}
$$
</div>

With more observed data ($n \rightarrow \infty$), the posterior mean is dominated by the sample mean $\bar{y}$.

The `posterior precision decomposition` is:

<div>
$$
    \tilde{\tau}_n^2 = \tilde{\tau}_0^2 + n\tilde{\sigma}^2
$$
</div>

This says the "posterior information" is equal to the prior information plus information from the data.

## Posterior predictive checks

A common way to assess model fit is to perform a `posterior predictive check`. The idea is our posterior for $\theta$ represents our updated beliefs about the mean carapace length after observing data. If we appropriately modeled carapace lengths (by our normal model), then **new predictions** from our sampling model should look similar to our observed data.

Recall what our observed data look like in Figure 1. When we generate new predictions, we want to see the new values to be within 10 to 50, roughly symmetric, and center around 30.

The new predictions are generated from our **sampling model**, because the sampling model describes how our data is generated given parameters $\theta$ and $\sigma^2$.

<div>
$$
    \tilde{Y} \mid \theta, \sigma^2 \sim N(\theta, \sigma^2 = 45)
$$
</div>

Note that the tilde on $Y$ represents a new prediction. Since the posterior mean $E(\theta \mid y)$ represents the mean updated belief about $\theta$, we could plug its value (32.1mm) in and generate predictions from:

<div>
$$
    \tilde{Y} \sim N(32.1, 45)
$$
</div>

However, this method doesn't account for the posterior uncertainty in $\theta$. A better way is to **average** over the posterior distribution of $\theta$ by computing the `posterior predictive distribution`:

<div>
$$
    p(\tilde{y} \mid y) = \int \underbrace{p(\tilde{y} \mid \theta, \sigma^2 = 45)}_{\text{prediction model}} \underbrace{p(\theta \mid y, \sigma^2 = 45)}_{\text{posterior for } \theta} d\theta
$$
</div>

Formally, the posterior predictive distribution can be difficult to compute in closed-form. However, it is much simple to sample from the posterior prediction distribution using Monte Carlo. The process is we first sample a $\theta$ from the posterior distribution, and then draw from the predictive model centered at $\theta$:

<div>
$$
    \begin{gathered}
        \text{draw } \theta^{(1)} \sim p(\theta \mid y, \sigma^2=45) \implies \text{draw } \tilde{y}^{(1)} \sim p(\tilde{y} \mid \theta^{(1)}, \sigma^2=45) \\
        \text{draw } \theta^{(2)} \sim p(\theta \mid y, \sigma^2=45) \implies \text{draw } \tilde{y}^{(2)} \sim p(\tilde{y} \mid \theta^{(2)}, \sigma^2=45) \\
        \vdots \\
        \text{draw } \theta^{(S)} \sim p(\theta \mid y, \sigma^2=45) \implies \text{draw } \tilde{y}^{(S)} \sim p(\tilde{y} \mid \theta^{(S)}, \sigma^2=45)
    \end{gathered}
$$
</div>

Then the set of values `$\{\tilde{y}^{(1)}, \cdots, \tilde{y}^{(S)}\}$` is the empirical distribution of the posterior predictive. For our crab example, the densities of the observed data and the posterior predictive are given below.

```r
# Prior and posterior parameters
sigma <- sqrt(45)
mu_0 <- 28
tau_0 <- 8
n <- nrow(crabs)
mu_n <- ((1/tau_0^2)*mu_0 + (n/sigma^2)*mean(crabs$CL)) /
  (1/tau_0^2 + n/sigma^2)
tau_n <- sqrt(1/((1/tau_0^2)+(n/sigma^2)))
S <- 1000

# MC samples from the posterior and posterior predictive distributions
set.seed(42)
theta_post <- rnorm(S, mean = mu_n, sd = tau_n)
theta_postpred <- rnorm(S, mean = theta_post, sd = sigma)
```

The two densities look similar and are centered around the same value. However, the observed data density has heavier tails (i.e. is more spread out) than the posterior predictive, indicating our sample model underestimated the variability $\sigma^2$. We might assume a larger value for $\sigma^2$ to better represent the variability in the data, or instead of fixing $\sigma^2$ we place a prior on both parameters and get a posterior distribution for both. This would allow the data to drive what the value should be.

{{< figure src="posterior_pred.png" caption="Comparison of the density of the observed data and the posterior predictive distribution." numbered="true" >}}

## Improper priors

One other issue that hasn't arisen yet is known as an `improper prior`. Suppose we have the normal sampling model with unknown mean $\theta$ and known variance $\sigma^2$. What if we have no prior information about $\theta$?

A candidate prior model is:

<div>
$$
    p(\theta) = 1, \quad -\infty < \theta < \infty
$$
</div>

where any real number is equally plausible. This prior model is _not_ a valid distribution mathematically. However, using this expression as the prior results in a closed-form posterior which _is_ a valid distribution!

Note that improper priors can result in model instability. For this particular model[^linear-regression], this prior works because we have a large amount of data. With a smaller dataset, suggesting such a large range of possible values of $\theta$ can make things messy. In general, we should be cautious when using them, and should use a valid prior distribution with large variance if possible.

[^linear-regression]: The other place we might see improper priors is in linear regression. Often times we place improper priors on the regression coefficients.
