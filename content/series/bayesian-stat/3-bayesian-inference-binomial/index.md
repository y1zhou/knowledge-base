---
title: "Bayesian Inference for the Binomial Model"
date: 2021-01-25T10:20:00-05:00
summary: "The general procedure for Bayesian analysis. We use two different prior models and compare the resulting posteriors (visually and mathematically)." # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-bayesian-inference-binomial"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 30 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Person working on a test." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "s9CC2SKySJM" # Unsplash ID of the picture
---

In this lecture we will talk about the same model as in last week, but from a Bayesian perspective. A lot of the details (calculations, choices of priors) are omitted but will be discussed in the next lecture.

Let $A$ be the event of interest from an experiment. Recall that the frequentist view of probability is the limiting relative frequency interpretation:

<div>
$$
P(A) = \lim_{n \rightarrow \infty} \frac{n(A)}{n}
$$
</div>

and the Bayesian takes the "belief" statement:

<div>
$$
P(A) = \text{degree of belief that event } A \text{ will occur.}
$$
</div>

As we will see, this difference in interpretations has ramifications on how statistical inference is conducted.

The example we use comes from _Hoff, Section 3.1_. A sample of 129 females aged 65 or older were asked whether or not they were generally happy in the 1998 General Social Survey. Of these, 118 reported being generally happy, while the remaining 11 reported not being generally happy. Some notations:

<div>
$$
\begin{gathered}
    n = 129 = \text{sample size} \\
    y = 118 = \text{# who reported being generally happy}
\end{gathered}
$$
</div>

Our goal is to estimate the true rate $\theta$ of generally happy females aged 65 or older.

## Ingredients for Bayesian modeling

Every Bayesian analysis requires specification of the following two pieces:

-   `Sampling model`: $p(y \mid \theta)$. This represents how likely is it that we observed $y$ females generally happy in our sample of size $n$, if the true happiness rate is $\theta$? This is going to be identical to how we treat it in the frequentist setting. We're relating our sample $y$ to the population parameter $\theta$.

-   `Prior model`: $p(\theta)$ summarizes beliefs associated with different values of the true happiness rate $\theta$, **prior** to observing data $y$. This doesn't make sense in the frequentist setting because the parameters would be fixed.

Both of these require us to specify distributions for these models. Once we have the sampling model and the prior model, our goal is to study the `posterior distribution` denoted $p(\theta \mid y)$. The posterior distribution summarizes beliefs associated with different values of the true happiness rate $\theta$, **after** observing data $y$. It's a procedure of updating our prior beliefs in a way that's consistent with the collected data $y$.

Given the <span style="color: crimson">prior</span> and the <span style="color: orange">sampling model</span>, we can compute the posterior by `Bayes' theorem`:

<div>
$$
p(\theta \mid y) = \frac{\textcolor{orange}{p(y \mid \theta)} \textcolor{crimson}{p(\theta)}}{p(y)}
$$
</div>

Here $p(y)$ is the `normalizing constant`, and we can compute this using just $p(y \mid \theta)$ and $p(\theta)$. In the next few lectures we will focus on posteriors that can be directly computed, but there are complications that arise from this expression which we will get to later.

### Sampling model

Let $Y$ be a random variable counting the number of "successes" (e.g. generally happy females) within a sample of size $n$. If $\theta$ represents the true rate of "success" in the population, then an appropriate sampling model is the binomial model:

<div>
$$
Y \mid \theta \sim Bin(n, \theta)
$$
</div>

The probability mass function is

<div>
$$
\begin{equation}\label{eq:binom-pdf}
    p(y \mid \theta) = \binom{n}{y} \theta^y (1-\theta)^{n-y}, \quad y = 0, 1, \cdots, n
\end{equation}
$$
</div>

Note this choice is identical to how we would traditionally specify a model, i.e. if we were doing frequentist inference.

### Non-informative prior model 1

Suppose we don't have any real prior beliefs about $\theta$, in other words the only prior belief we have about $\theta$ is that any value between 0 and 1 is equally plausible. An appropriate, `non-informative` prior model is:

<div>
$$
\theta \sim Unif(0, 1)
$$
</div>

The probability density function of our uniform distribution is

<div>
$$
p(\theta) = 1, \quad 0 < \theta < 1
$$
</div>

Again, because we are now interpreting probabilities as degrees of belief, it makes sense (and is necessary) to talk about probability distributions for $\theta$, whereas this would have no merit from a frequentist perspective.

Now we can use Bayes' Theorem to directly compute the posterior distribution. Note that the normalizing constant is

<div>
$$
p(y) = \int_{-\infty}^\infty p(y \mid \theta) p(\theta) d\theta
$$
</div>

Using the `kernel trick`[^kernel-trick], we will show that the posterior is a beta distribution:

[^kernel-trick]: After writing out the posterior distribution using Bayes' Theorem, we make the denominator look like the integrals of a pdf of some known distribution so that we don't have to solve it.

<div>
$$
\left.\begin{aligned}
    y \mid \theta \sim Bin(n, \theta) \\
    \theta \sim Unif(0, 1)
\end{aligned}\right\} \Rightarrow \theta \mid y \sim Beta(y+1, n-y+1)
$$
</div>

Here we skip the calculations because the uniform prior is a special case of the later beta prior. The beta distribution is a continuous distribution with two parameters, which is defined on the interval from 0 to 1. For our specific data, we have $n = 129$ and $y = 118$, so

<div>
$$
    \theta \mid_{y=118} \sim Beta(118+1, 129-118+1) = Beta(119, 12)
$$
</div>

#### Visualizing the model

We can visualize how the observed data updates prior 1 beliefs into posterior 1 beliefs about $\theta$[^report-bayesian]:

[^report-bayesian]: When we report Bayesian results, it's always good to plot the prior and posterior simultaneously because this gives us all the information about our beliefs about $\theta$ before and after observing the data.

```r
library(dplyr)
library(tidyr)
library(ggpubr)

# Female happiness data
n <- 129
y <- 118

# Uniform prior for theta (yielding posterior 1)
# Evaluate prior and posterior pdfs on grid
theta <- seq(0, 1, length.out = 1000)
prior1_pdf <- dunif(grid, 0, 1)
post1_pdf <- dbeta(grid, y + 1, n - y + 1)

# Plot the pdf
data.frame(theta, Prior = prior1_pdf, Posterior = post1_pdf) %>%
  pivot_longer(
    -theta,
    names_to = "Distribution",
    values_to = "Density"
  ) %>%
  mutate(Distribution = as.factor(Distribution),
         Distribution = relevel(Distribution, "Prior")) %>%
  ggline(x = "theta", y = "Density", color = "Distribution",
         xlab = expression(theta), ylab = "Density",
         palette = "npg", plot_type = "l", numeric.x.axis = T)
```

{{< figure src="noninformative_prior_pdf.png" caption="Non-informative prior 1 and posterior 1 for happiness rate $\theta$." numbered="true" >}}

We can see the distribution looks very different. One way to interpret this is our expected true rate of happiness is about 90% based off the observed data. We can also come up with suitable range of values for our guess of $\theta$. We attribute very low degrees of beliefs to values below 80%.

### Informative prior model 2

If we think about the previous prior, it might not make sense to say that "only 1% of females are generally happy". The statement seems implausible even before observing any data, so we might want to specify a different prior distribution.

What if we had knowledge of results from a previous General Social Survey, which identified that 45 of 60 females aged 65 or older were generally happy? An `informative prior` based on this information may satisfy the following:

-   Center: based on this previous study, we might guess $\theta$ is relatively close to $45 / 60 = 0.75$.
-   Spread: if 75% is plausible, how much "less plausible" is 76%, i.e. how quickly does the prior distribution go down towards zero? This might depend on the sample size of the previous study -- a larger sample size means we might place more trust in the prior beliefs, and that results in less spread in the prior distribution.

An appropriate choice of prior is $\theta \sim Beta(a = 45, b = 15)$ where $a$ and $b$ are `hyperparameters`. We can choose different hyperparameters[^hyperparam-interpretation] to adapt to our prior beliefs about $\theta$. In this case, our prior distribution is centered around 75%.

[^hyperparam-interpretation]: In this specific case, the prior hyperparameters have nice interpretations: $a$ is the prior number of successes, $b$ is the prior number of failures, and $a+b$ is the prior sample size.

We will show that this choice of sampling and prior model once again results in a beta posterior. Recall from $\eqref{eq:binom-pdf}$ that our sampling distribution $p(y \mid \theta)$ is binomial. The pdf of prior 2 is

<div>
$$
    \begin{equation}\label{eq:beta-pdf}
        p(\theta) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \theta^{a-1} (1-\theta)^{b-1}, \quad 0 < \theta < 1
    \end{equation}
$$
</div>

where

<div>
$$
    \Gamma(z) = \int_0^\infty x^{z-1} e^{-x} dx
$$
</div>

is the `gamma function`. Using Bayes' Theorem, the posterior distribution is

<div>
$$
    \begin{align}
        p(\theta \mid y) &= \frac{p(y \mid \theta) p(\theta)}{p(y)} = \frac{\textcolor{crimson}{p(y \mid \theta)} \textcolor{orange}{p(\theta)}}{\int_{-\infty}^\infty \textcolor{crimson}{p(y \mid \theta)} \textcolor{orange}{p(\theta)} d\theta} \\
        &= \cfrac{\textcolor{crimson}{\binom{n}{y}\theta^y (1-\theta)^{n-y}} \textcolor{orange}{\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \theta^{a-1} (1-\theta)^{b-1}}}{\int_{-\infty}^\infty \textcolor{crimson}{\binom{n}{y}\theta^y (1-\theta)^{n-y}} \textcolor{orange}{\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)} \theta^{a-1} (1-\theta)^{b-1}} d\theta} \label{eq:posterior2}
    \end{align}
$$
</div>

If we look at the denominator, some constant terms that don't depend on $\theta$ can be taken out of the integral:

<div>
$$
    \begin{align}
        \cdots &= \textcolor{green}{\binom{n}{y}\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}} \int_{-\infty}^\infty \theta^y (1-\theta)^{n-y} \theta^{a-1} (1-\theta)^{b-1} d\theta \\
        &= \textcolor{green}{c} \int_{-\infty}^\infty \theta^{y+a-1} (1-\theta)^{n-y+b-1} d\theta \label{eq:posterior2-denom}
    \end{align}
$$
</div>

Now let `$\tilde{a} = y+a$` and `$\tilde{b} = n-y+b$`, we have

<div>
$$
    \begin{align}
        \eqref{eq:posterior2-denom} &= c \int_{-\infty}^\infty \theta^{\tilde{a}-1} (1-\theta)^{\tilde{b}-1} d\theta \\
        &= c \cdot \frac{\Gamma(\tilde{a})\Gamma(\tilde{b})}{\Gamma(\tilde{a}+\tilde{b})} \cdot \frac{\Gamma(\tilde{a}+\tilde{b})}{\Gamma(\tilde{a})\Gamma(\tilde{b})} \int_{-\infty}^\infty \theta^{\tilde{a}-1} (1-\theta)^{\tilde{b}-1} d\theta \\
        &= c \cdot \frac{\Gamma(\tilde{a})\Gamma(\tilde{b})}{\Gamma(\tilde{a}+\tilde{b})} \int_{-\infty}^\infty \frac{\Gamma(\tilde{a}+\tilde{b})}{\Gamma(\tilde{a})\Gamma(\tilde{b})} \theta^{\tilde{a}-1} (1-\theta)^{\tilde{b}-1} d\theta \label{eq:posterior2-denom-final}
    \end{align}
$$
</div>

See that the integrand in `$\eqref{eq:posterior2-denom-final}$` takes the same form as the pdf of the beta distribution in `$\eqref{eq:beta-pdf}$`. Thus the entire integral part evaluates to 1, and the denominator is just a constant term. Going back to the posterior model,

<div>
$$
    \begin{aligned}
        \eqref{eq:posterior2} &= \cfrac{c \theta^{y+a-1} (1-\theta)^{n-y+b-1}}{c \cdot \frac{\Gamma(\tilde{a})\Gamma(\tilde{b})}{\Gamma(\tilde{a}+\tilde{b})}} \\
        &= \frac{\Gamma(\tilde{a}+\tilde{b})}{\Gamma(\tilde{a})\Gamma(\tilde{b})} \theta^{\tilde{a}-1} (1-\theta)^{\tilde{b}-1} \\
        &= \text{pdf of } Beta(\tilde{a}, \tilde{b})
    \end{aligned}
$$
</div>

So with the informative prior, we have

<div>
$$
    \left.\begin{aligned}
        y \mid \theta \sim Bin(n, \theta) \\
        \theta \sim Beta(a, b)
    \end{aligned} \right\} \implies \theta \mid y \sim Beta(y+a, n-y+b)
$$
</div>

We can also see that the uniform prior we originally specified is equivalent to a $Beta(a=1, b=1)$ distribution, thus this result is a generalization of what we stated for posterior 1.

#### Visualizing the model

With a slight modification of the code above, we can visualize how the observed data updates our informative prior 2 beliefs into our posterior 2 beliefs about $\theta$:

{{< figure src="informative_prior_pdf.png" caption="Informative prior 2 and posterior 2 for happiness rate $\theta$." numbered="true" >}}

### Shortcut for computing posterior

We compute the posterior distribution via Bayes' Theorem. Computing the normalizing constant (denominator) is a key part of deriving the posterior. We used the [kernel](<https://en.wikipedia.org/wiki/Kernel_(statistics)#Bayesian_statistics>) trick to evaluate the integral, by recognizing that the integrand looked like a $Beta(y+a, n-y+b)$ density function.

In many cases, the normalizing constant is difficult, if not impossible, to compute directly. But do we actually _need to_ compute this constant? Suppose we are interested in estimating $\theta$ by either `$\theta_1=0.3$` or `$\theta_2 = 0.9$`. Given data $y$, we want to know which value is a better estimate. One way to think about this is we would prefer `$\theta_1$` if

<div>
$$
    p(\theta_1 \mid y) > p(\theta_2 \mid y) \implies \frac{p(\theta_1 \mid y)}{p(\theta_2 \mid y)} > 1
$$
</div>

This is like a likelihood ratio, but in terms of posterior densities. Note that

<div>
$$
    \frac{p(\theta_1 \mid y)}{p(\theta_2 \mid y)} = \frac{p(y \mid \theta_1) p(\theta_1) / p(y)}{p(y \mid \theta_2) p(\theta_2) / p(y)} = \frac{p(y \mid \theta_1) p(\theta_1)}{p(y \mid \theta_2) p(\theta_2)}
$$
</div>

does **not** dependent on the normalizing constant $p(y)$. This means {{<hl>}}the posterior's shape only depends on the numerator.{{</hl>}} If we plot the true posterior and the numerator side by side, they will be exactly the same shape. $p(y)$ only adjusts how tall the peak is so that the density integrals to 1. We'd prefer `$\theta_2$` over `$\theta_1$` regardless of which plot we look at.

{{< figure src="numerator_shape.png" caption="The shape of the posterior is the same as the numerator." numbered="true" >}}

In light of this, a simpler way to compute the posterior is to work with `proportionality statements`:

<div>
$$
    \begin{aligned}
        p(\theta \mid y) &= \frac{p(y \mid \theta) p(\theta)}{\textcolor{crimson}{p(y)}} \propto p(y \mid \theta) p(\theta) \\
        &\propto \textcolor{crimson}{\binom{n}{y}} \theta^y (1-\theta)^{n-y} \textcolor{crimson}{\frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)}} \theta^{a-1} (1-\theta)^{b-1} \\
        &\propto \theta^{y+a-1} (1-\theta)^{n-y+b-1} \implies \text{kernel of } Beta(y+a, n-y+b)
    \end{aligned}
$$
</div>

All the colored terms are constant with respect to $\theta$. Thus, the posterior _is_ a $Beta(y+a, n-y+b)$ distribution. We will use this strategy throughout the rest of the course, where instead of evaluating the denominator, we will just focus on proportionality statements with the numerator and see if it looks like a distribution we recognize.

## Conjugacy

For the binomial model with a beta prior, we were able to *directly* recognize the posterior distribution. Sadly this doesn't always happen. Suppose we specified the following prior model instead:

<div>
$$
    p(\theta) = \frac{1}{\sqrt{2\pi\sigma^2}} e^{-\frac{(\theta-\mu)^2}{2\sigma^2}}, \quad -\infty < \theta < \infty
$$
</div>

for some small $\sigma$. This corresponds to the normal distribution that is centered at $\mu$, and $\sigma$ controls the amount of spread. Now if we try to derive the posterior distribution,

<div>
$$
    \begin{aligned}
        p(\theta \mid y) & \propto p(y \mid \theta) p(\theta) \\
        &\propto \textcolor{crimson}{\binom{n}{y}} \theta^y (1-\theta)^{n-y} \textcolor{crimson}{\frac{1}{\sqrt{2\pi\sigma^2}}} e^{-\frac{(\theta-\mu)^2}{2\sigma^2}} \\
        &\propto \theta^y (1-\theta)^{n-y} e^{-\frac{(\theta-\mu)^2}{2\sigma^2}}
    \end{aligned}
$$
</div>

This is no longer the kernel of some distribution that we recognize. This doesn't mean our prior belief or the sampling model is wrong. We'd still get a valid posterior distribution, but we won't be able to find the normalizing constant in closed form. In later lectures we'll talk more about how to handle situations like this.

The prior $\theta \sim Beta(a, b)$ was conveniently chosen because we obtain a closed-form posterior distribution. This choice of convenience is known as a conjugate prior. Formally speaking, let $\mathcal{F}$ be the class of sampling models $p(y \mid \theta)$ and $\mathcal{P}$ be the class of prior distributions $p(\theta)$. We say the prior class $\mathcal{P}$ is `conjugate` for sampling model class $\mathcal{F}$ if:

<div>
$$
    p(\theta) \in \mathcal{P} \text{ and } p(y \mid \theta) \in \mathcal{F} \implies p(\theta \mid y) \in \mathcal{P}.
$$
</div>

Intuitively, the prior $p(\theta)$ and posterior $p(\theta \mid y)$ come from the same family of distributions. For example, in the binomial model we have

<div>
$$
    \begin{aligned}
        &\text{Prior}: \theta \sim Beta(a, b) \\
        &\text{Posterior}: \theta \mid y \sim Beta(y+a, n-y+b)
    \end{aligned}
$$
</div>

We can say the beta prior is conjugate for the binomial model, since our posterior is also beta-distributed.

### Pros and Cons

Pros to conjugate prior specification:

- Allows direct computation of posterior, i.e. a closed-form, recognizable posterior distribution.
- Easy to understand and very interpretable.
- Good starting choice.

Cons to conjugate prior specification:

- Now always available for every model, particularly if $\theta$ consists of a large number of parameters (high-dimensional).
- Restrictive in terms of modeling choices. It might not always be the choice which best reflects prior beliefs about $\theta$.

## Numerical posterior summaries

It is useful to produce numerical summaries of the posterior distribution:

<div>
$$
    \theta \mid y \sim Beta(y+a, n-y+b)
$$
</div>

### Measures of center

First, we focus on measures of center. The `posterior mean` tells us what is the expected proportion of generally happy females aged 65 or older, given the observed data.

<div>
$$
    E(\theta \mid y) = \frac{a+y}{a+b+n} = \boxed{\frac{a+b}{a+b+n}} \cdot \textcolor{crimson}{\frac{a}{a + b}} + \boxed{\frac{n}{a+b+n}} \cdot \textcolor{orange}{\frac{y}{n}}
$$
</div>

The equation above can be interpreted as a weighted sum of the <span style="color: crimson">prior mean for $\theta$</span> and the <span style="color: orange">observed sample proportion</span>. The boxed terms are weights on the prior mean and the observed data, respectively.

The posterior mean is a weighted sum of the prior mean and the sample mean. As the sample size `$n \rightarrow \infty$`, `$E(\theta \mid y) \rightarrow y/n$`, i.e. the data "overrides" the prior information.

An alternative measure of center is the `maximum a posteriori (MAP)` estimate. This is also known as the `posterior mode`, and tells us what is the most likely proportion of generally happy females aged 65 or older, given the observed data[^map-estimate].

[^map-estimate]: This value corresponds to the peak in the plot of the posterior density function.

<div>
$$
    \hat\theta_{MAP} \triangleq \text{value of } \theta \text{ which maximizes the posterior } p(\theta \mid y)
$$
</div>

For our posterior,

<div>
$$
    \hat\theta_{MAP} = \frac{y+n-1}{n+a+b-2}
$$
</div>

We might report the MAP estimate over the posterior mean because the MAP is the value of $\theta$ for which our posterior beliefs are highest. This is especially useful for skewed distributions.

### Measures of variation

We may also be interested in the amount of variation in the posterior of $\theta$. The `posterior variance` measures how much uncertainty there is in the proportion of generally happy females aged 65 or older, given the observed data.

<div>
$$
    Var(\theta \mid y) = \frac{(y+a)(n-y+b)}{(n+a+b)^2 (n+a+b+1)}
$$
</div>

As the observed sample size $n$ increases, the posterior variation decreases. This measure however does not have the nice interpretation the posterior mean has. A better way to quantify uncertainty in our posterior beliefs about $\theta$ is by computing a _range of values_ for which it is likely $\theta$ falls within.

A `$100(1-\alpha)\%$` `posterior credible interval` for $\theta$ is determined by endpoints $l(Y)$ and $u(Y)$ which satisfy

<div>
$$
    \begin{equation}\label{eq:credible-interval}
        P\big( l(Y) < \theta < u(Y) \mid Y=y \big) = 1 - \alpha
    \end{equation}
$$
</div>

This is saying the posterior probability that $\theta$ is in the interval `$[l(Y), u(Y)]$` is equal to $1-\alpha$. This is actually a really nice interpretation because all we need to do is to find the two endpoints on the posterior distribution which we directly have.

Note that this is **not** the same as a confidence interval! A `$100(1-\alpha)\%$` confidence interval for $\theta$ would satisfy

<div>
$$
    P\big( l(Y) < \theta < u(Y) \big) = 1 - \alpha
$$
</div>

The frequentist CI works regardless of what sample we took from the population. More specifically, it is not conditional on the specific observed data because it's based off repeatedly generated data from the population. That being said, credible intervals and confidence intervals often coincide with each other. In the Bayesian setting, we always report credible intervals.

### Types of credible intervals

Now the question is how do we select $l(Y)$ and $u(Y)$. There are two different ways to construct intervals by determining these endpoints: a `percentile-based interval`, or a `highest posterior density (HPD) region`.

#### Percentile-based interval

{{< figure src="percentile_interval.png" caption="Zoomed-in posterior 1 for happiness rate $\theta$ and the 95% credible interval." numbered="true" >}}

If we choose $l(Y)$ and $u(Y)$ as percentiles of the posterior, then we will have $\eqref{eq:credible-interval}$. We want $1-\alpha$ in the shaded region. Choosing percentiles means we want `$\frac{\alpha}{2}$` probabilities in the left and right tails. The cutoff points can be found by

<div>
$$
    \begin{gathered}
        l(Y) = 100 \left(\frac{\alpha}{2}\right)^{th} \text{ percentile of } p(\theta \mid y) = 0.854 \\
        u(Y) = 100 \left(1 - \frac{\alpha}{2}\right)^{th} \text{ percentile of } p(\theta \mid y) = 0.951
    \end{gathered}
$$
</div>

The numbers are found in R using the `qbeta()` function. For example, the 2.5-th percentile $l(Y)$ can be found by `qbeta(0.025, y+1, n-y+1)`. We could say that with probability 0.95, the true rate $\theta$ is between 0.854 and 0.951. This is a much more straightforward interpretation than the frequentist interpretation of confidence intervals.

#### HPD region

The percentile-based interval is easy to compute and is most commonly used. However, some values of $\theta$ outside of a percentile-based CI might have higher posterior densities than the values within (due to the "asymmetry").

The HPD region is constructed in a way to resolve this issue. $l(Y)$ and $u(Y)$ are chosen to satisfy:

1. $P(l(Y) < \theta < u(Y) \mid Y=y) = 1-\alpha$.
2. Any value of $\theta$ between $l(Y)$ and $u(Y)$ must have higher posterior density than values outside the interval.

In R, the `HDInterval::hdi()` function can be used to find the HPD region. As a comparison to the 95% percentile CI of (0.854, 0.951), the HPD region is (0.858, 0.955). In this case the two are not very different because the posterior is quite symmetric, but the HPD region comes handy when the posterior is skewed or multimodal.

## Comparing the two posteriors

Now we can summarize and compare posteriors 1 and 2:

| Posterior          | `$\theta \mid y \sim Beta(119, 12)$` | `$\theta \mid y \sim Beta(163, 26)$` |
| :----------------- | -----------------------------------: | -----------------------------------: |
| Posterior mean     |                                0.908 |                                0.862 |
| Posterior variance |                               0.0006 |                               0.0006 |
| MAP                |                                0.915 |                                0.866 |
| 95% posterior CI   |                       (0.854, 0.951) |                       (0.810, 0.908) |

Posterior 2 is shifted to the left of 1, since we specified an informative prior which was centered around $\theta = 0.75$. Both posteriors are concentrated on a small set of values.

## Sensitivity analysis

Re-visiting our prior specification, suppose we believed from past studies that the true rate of generally happy females aged 65 or older was closer to 75%. Consider the conjugate prior for $\theta$ ($\theta \sim Beta(a, b)$), with the following hyperparameter choices:

| Prior | "Prior" sample size |  $a$ |  $b$ |
| :---- | ------------------: | ---: | ---: |
| 1     |                   2 |    1 |    1 |
| 2     |                   4 |    3 |    1 |
| 3     |                  60 |   45 |   15 |
| 4     |                 120 |   90 |   30 |

Recall that the prior mean of $\theta$ is $E(\theta) = a/(a+b)$. We can interpret $a$ as the "prior" number of successes, and $b$ as the "prior" number of failures. We can control $a$ and $b$ so that the prior mean is centered around 0.75. If we increase the prior sample size, we put more emphasize on the prior beliefs.

In the table above, prior 1 is the non-informative prior, and the other three are centered at $\theta = 0.75$. The posterior distribution under these prior choices are:

<div>
$$
    \begin{gathered}
        (1): \theta \mid y=118 \sim Beta(119, 12) \\
        (2): \theta \mid y=118 \sim Beta(121, 12) \\
        (3): \theta \mid y=118 \sim Beta(163, 26) \\
        (4): \theta \mid y=118 \sim Beta(206, 41)
    \end{gathered}
$$
</div>

A study of how sensitive the posterior is to our prior assumptions is called a `sensitivity analysis`. We're looking at how sensitive is our posterior distribution to different reasonable choices of prior distributions, and it is remarkably important for *any* Bayesian analysis.

As we go from (1) to (4), we place more emphasis on our prior beliefs. The 95% HPD credible intervals for $\theta$ are:

| Prior | Credible interval |
| ----- | ----------------- |
| (1)   | (0.858, 0.955)    |
| (2)   | (0.860, 0.955)    |
| (3)   | (0.813, 0.910)    |
| (4)   | (0.789, 0.880)    |

The first two are very similar to each other due to their small sample sizes. As the sample size increases, the CIs shift to the left. The take-away message is if one does not have any prior information, using a non-informative prior to "let the data drive inference" is a good choice. If one does have prior information, use it carefully especially if the posterior is sensitive to the choice of prior, i.e. if there's not much observed data.
## Concluding remarks

1. Every Bayesian analysis proceeds in a similar fashion to this one. In the coming lectures, we will talk more about different details associated with this process.
2. Neither posterior nor 2 should be deemed correct or incorrect. It all depends on the amount of prior information available. In this example, the posterior was quite different between the two prior choices. It is often beneficial to study how sensitive the posterior is to prior information.
3. Generally, the prior distribution should **not** depend on the observed data. This is often referred to as "double-dipping", and can lead to faulty conclusions due to over-reliance on the data. You can base the prior off previous studies, scientific knowledge, or other intuition.
