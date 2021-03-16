---
title: "The Normal Model in a Two Parameter Setting"
date: 2021-03-15T10:20:00-05:00
summary: "This lecture discusses Bayesian inference of the normal model, particularly the case where we are interested in joint posterior inference of the mean and variance simultaneously. We discuss approaches to prior specification, and introduce the Gibbs sampler as a way to generate posterior samples if full conditional distributions of the parameters are available in closed-form." # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-normal-two-param-setting"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 70 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

Starting from this lecture, we shift gears from looking at posteriors in a mathematical perspective to the computational aspect of Bayesian inference. We are going to revisit the normal model, but this time in a two parameter setting. This is a natural model to consider for the first computational method we talk about -- Gibbs sampling.

## Motivating example

So far, we have focused on inference of a single parameter $\theta$. In addition, we have primarily specified conjugate prior distributions so that we get a closed-form posterior for $\theta$:

| Sampling model | Parameter                      | Prior         | Posterior     |
| -------------- | ------------------------------ | ------------- | ------------- |
| Binomial       | $\theta$ (success probability) | Beta          | Beta          |
| Poisson        | $\theta$ (mean)                | Gamma         | Gamma         |
| Normal         | $\theta$   (mean)              | Normal        | Normal        |
| Normal         | $\sigma^2$ (variance)          | Inverse-Gamma | Inverse-Gamma |

For the normal model specifically, we were either estimating the mean given that the variance was fixed, or estimating the variance given that the mean was fixed. A natural question is what if we are interested in inference of *multiple parameters*? We want to assume that both $\theta$ and $\sigma^2$ are unknown.

Suppose we are interested in estimating the mean and variance in calories per serving of cereal. The data set with $n = 77$ cereals come from [Kaggle](https://www.kaggle.com/crawford/80-cereals). The data has a mean of $\bar{y} = 106.88$ and a standard deviation of $s = 19.48$.

{{< figure src="calories_hist.png" caption="Histogram of calories per serving of the cereals." numbered="true" >}}

The data looks roughly symmetric and is centered around 105. It has somewhat heavier tails than a typical normal distribution.

## Sampling model

Let $Y_i$ be a random variable representing the calories per serving for the $i$-th cereal for $i = 1, \cdots, n$ where $n=77$. Then, assuming conditional independence, an appropriate sampling model is:

<div>
$$
    Y_1, \cdots, Y_n \mid \theta, \sigma^2 \overset{i.i.d.}{\sim} N(\theta, \sigma^2)
$$
</div>

The parameters of interest are:

- $\theta$, the true mean cereal calories per serving, and
- $\sigma^2$, the true variance of cereal calories per serving.

The joint probability density function is given by:

<div>
$$
    p(y \mid \theta, \sigma^2) = \prod_{i=1}^n p(y_i \mid \theta, \sigma^2) = (2\pi\sigma^2)^{-\frac{n}{2}} \exp \left\{ -\frac{1}{2\sigma^2} \sum_{i=1}^n (y_i - \theta)^2 \right\}
$$
</div>

Nothing new so far!

## Prior model

We want to estimate two parameters, $\theta$ and $\sigma^2$, *simultaneously*. In earlier lectures we assumed one of the parameters was fixed and known, and only estimated the other. When $\theta$ is unknown and $\sigma^2$ is known, we had:

<div>
$$
    \begin{gathered}
        \text{Prior: } \theta \sim N(\mu_0, \tau_0^2) \\
        \text{Posterior: } \theta \mid y \sim N(\mu_n, \tau_n^2)
    \end{gathered}
$$
</div>

When $\sigma^2$ is unknown and $\theta$ is known:

<div>
$$
    \begin{gathered}
        \text{Prior: } \sigma^2 \sim IG(a, b) \\
        \text{Posterior: } \sigma^2 \mid y \sim IG(a_n, b_n)
    \end{gathered}
$$
</div>

In our case, both parameters are unknown now, and we must specify a prior distribution for the parameter set $(\theta, \sigma^2)$.

### Choice 1

Our first choice is to use a **conditional specification**. The multiplication rule in probability shows that:

<div>
$$
    P(A, B) = P(A \mid B) P(B)
$$
</div>

So similarly we can specify the prior distribution for the two parameters conditionally:

<div>
$$
    p(\theta, \sigma^2) = p(\theta \mid \sigma^2) p(\sigma^2),
$$
</div>

where

<div>
$$
    \sigma^2 \sim IG(a, b), \quad \theta \mid \sigma^2 \sim N(\mu_0, g_0(\sigma^2))
$$
</div>

Remember anytime we condition on something, we're treating the value as fixed. Here we're treating $\sigma^2$ as fixed, and the variance of $\theta \mid \sigma^2$ depends on the value of $\sigma^2$. Specifying the prior in this way actually leads to `conditional conjugacy` in the sense that we would get the posterior in the form of:

<div>
$$
    \theta \mid \sigma^2, y \sim N(\mu_n, g_n(\sigma^2))
$$
</div>

We will not discuss this case in detail as it does not seem to be used much in practice.

### Choice 2

Our focus is on the **independent specification**. Assuming *a priori* that $\theta$ and $\sigma^2$ are independent of each other,

<div>
$$
    p(\theta, \sigma^2) = p(\theta) p(\sigma^2),
$$
</div>

then we could specify the same priors we discussed before:

<div>
$$
    \begin{gathered}
        \theta \sim N(\mu_0, \tau_0^2) \\
        \sigma^2 \sim IG(a, b)
    \end{gathered}
$$
</div>

This is not constraining us to say $\theta$ must depend on $\sigma^2$ in the prior distribution, or vice versa. The hyperparameters can also be chosen in the same ways as previously discussed.

## Computing the joint posterior

The `joint prior probability density function` is:

<div>
$$
    \begin{aligned}
        p(\theta, \sigma^2) = p(\theta) p(\sigma^2) &= \left[ (2\pi\tau_0^2)^{-\frac{1}{2}} \exp \left\{ -\frac{1}{2\tau_0^2} (\theta - \mu_0)^2 \right\} \right] \\
            &\quad \times \left[  \frac{b^n}{\Gamma(a)} \left( \frac{1}{\sigma^2} \right)^{a+1} \exp \left\{ -\frac{b}{\sigma^2} \right\} \right] \quad \text{for } \sigma^2 > 0
    \end{aligned}
$$
</div>

Under this independent prior specification, we can try to compute the joint posterior distribution:

<div>
$$
    \begin{aligned}
        p(\theta, \sigma^2 \mid y) &\propto p(y \mid \theta, \sigma^2) p(\theta, \sigma^2) \\
        &\propto (2\pi\sigma^2)^{-\frac{n}{2}} \exp \left\{ -\frac{1}{2\sigma^2} \sum_{i=1}^n (y_i - \theta)^2 \right\} \\
        &\quad \times \left[ (2\pi\tau_0^2)^{-\frac{1}{2}} \exp \left\{ -\frac{1}{2\tau_0^2} (\theta - \mu_0)^2 \right\} \right] \left[  \frac{b^n}{\Gamma(a)} \left( \frac{1}{\sigma^2} \right)^{a+1} \exp \left\{ -\frac{b}{\sigma^2} \right\} \right]
    \end{aligned}
$$
</div>

This is a huge mess and seems hopeless to simplify. We may exploit the multiplication rule for probabilities:

<div>
$$
    p(\theta, \sigma^2 \mid y) = p(\theta \mid \sigma^2, y)p(\sigma^2 \mid y)
$$
</div>

If we can obtain closed-form distributions for $p(\theta \mid \sigma^2, y)$ and $p(\sigma^2 \mid y)$ on the RHS, we could generate samples from the joint posterior conditionally via Monte Carlo sampling. Specifically,

<div>
$$
    \begin{gathered}
        \text{draw } \sigma^{2^{(1)}} \sim p(\sigma^2 \mid y) \implies \text{draw } \theta^{(1)} \sim p(\theta \mid \sigma^{2^{(1)}}, y) \\
        \text{draw } \sigma^{2^{(2)}} \sim p(\sigma^2 \mid y) \implies \text{draw } \theta^{(2)} \sim p(\theta \mid \sigma^{2^{(2)}}, y) \\
        \vdots \\
        \text{draw } \sigma^{2^{(S)}} \sim p(\sigma^2 \mid y) \implies \text{draw } \theta^{(S)} \sim p(\theta \mid \sigma^{2^{(S)}}, y)
    \end{gathered}
$$
</div>

Then `$\{ (\theta^{(1)}, \sigma^{2^{(1)}}), \cdots, (\theta^{(S)}, \sigma^{2^{(S)}}) \}$` are samples from the joint posterior distribution.

### Computations
The first piece we need is the posterior for $\theta$ conditional on $\sigma^2$:

<div>
$$
    \begin{aligned}
        p(\theta \mid \sigma^2, y) &\propto p(y \mid \theta, \sigma^2) p(\theta, \sigma^2) \\
        &\propto p(y \mid \theta, \sigma^2) p(\theta) p(\sigma^2)
    \end{aligned}
$$
</div>

Since we want the posterior for $\theta$, $p(\sigma^2)$ can be removed as it doesn't depend on $\theta$. The sampling model $p(y \mid \theta, \sigma^2)$ and the prior distribution $p(\theta)$ are both normal, so we would get

<div>
$$
    \theta \mid \sigma^2, y \sim N(\mu_n, \tau_n^2)
$$
</div>

This is exactly the same as the posterior in the unknown $\theta$, known $\sigma^2$ case, because conditioning on $\sigma^2$ means we assume its value is *known and fixed*.

Next we derive the posterior for $\sigma^2$:

<div>
$$
    p(\sigma^2 \mid y) \propto p(y \mid \sigma^2) p(\sigma^2)
$$
</div>

As our sampling model depended on $\theta$ **and** $\sigma^2$, we need to re-introduce $\theta$ into $p(y \mid \sigma^2)$:

<div>
$$
    \begin{aligned}
        p(\sigma^2 \mid y) &\propto p(y \mid \sigma^2) p(\sigma^2) \\
        &\propto \int p(y, \theta \mid \sigma^2)d\theta \cdot p(\sigma^2) \\
        &\propto \int p(y \mid \theta, \sigma^2) p(\theta) d\theta \cdot p(\sigma^2)
    \end{aligned}
$$
</div>

The integral is problematic for obtaining a closed-form distribution. As a consequence, we do not have a built-in R function which can directly generate Monte Carlo samples from $p(\sigma^2 \mid y)$.

### Full conditional distributions

While we can't obtain $p(\sigma^2 \mid y)$ in closed form, we **can** get a closed form distribution for $p(\sigma^2 \mid \theta, y)$:

<div>
$$
    \begin{aligned}
        p(\sigma^2 \mid \theta, y) &\propto p(y \mid \theta, \sigma^2) p(\theta, \sigma^2) \\
        &\propto p(y \mid \theta, \sigma^2) p(\theta) p(\sigma^2) \\
        &\propto \underbrace{p(y \mid \theta, \sigma^2)}_{\text{normal}} \underbrace{p(\sigma^2)}_{IG} \\
        \implies \sigma^2 \mid \theta, y &\sim IG(a_n, b_n)
    \end{aligned}
$$
</div>

Again, this is identical to the case where $\theta$ is fixed and $\sigma^2$ is unknown. This suggests that if we have data $y$ and we know the value of $\theta$, then we can generate a value for $\sigma^2$.

The distributions $p(\theta \mid \sigma^2, y)$ and $p(\sigma^2 \mid \theta, y)$ are called `full conditional distributions`, since they are the distribution of one parameter **given** everything else (i.e. the data and the other parameter) is fixed[^full-cond-dist]. The product of full conditional distributions is **not** the joint posterior we desire:

[^full-cond-dist]: Note that $p(\sigma^2 \mid y)$ is not a full conditional distribution, because it's not dependent on the value of $\theta$.

<div>
$$
    p(\theta, \sigma^2 \mid y) \neq \underbrace{p(\theta \mid \sigma^2, y)}_{(1)} \underbrace{p(\sigma^2 \mid \theta, y)}_{(2)}
$$
</div>

So why is this useful? Going back to Monte Carlo sampling and using the RHS of the expression above, if we know $\sigma^2$, then we can sample $\theta$ using (1); if we know $\theta$, then we can sample $\sigma^2$ using (2). We can still generate a **dependent** sequence of $(\theta, \sigma^2)$ values.

## Gibbs sampler

Suppose we have closed forms for the full conditional distributions of each unknown parameter, i.e. $p(\theta \mid \sigma^2, y)$ and $p(\sigma^2 \mid \theta, y)$. The following Monte Carlo sampling strategy is known as the `Gibbs sampler`:

1. Specify *starting values* for the parameters: $(\theta^{(1)}, \sigma^{2^{(1)}})$.
2. Suppose that at step $s$, the current values of the parameters are $(\theta^{(s)}, \sigma^{2^{(s)}})$. To generate values at step $s+1$:
   - Sample `$\theta^{(s+1)} \sim p(\theta, \sigma^{2^{(s)}}, y)$`, that is given the *current* value of $\sigma^2$.
   - Sample $\sigma^{2^{(s+1)}} \sim p(\sigma^2 \mid \theta^{(s+1)}, y)$. Here we sample the new value of $\sigma^2$ given the updated value of $\theta$.
3. Repeat step 2 $S$ times to generate a sequence of parameter values.

Unlike out previous methods of Monte Carlo sampling, parameter values at step $s+1$ depends on the previous step's values. This is known as the `Markov property`, where step $s+1$ only depends on the previous step $s$.

A question is does the sequence `$\{ (\theta^{(1)}, \sigma^{2^{(1)}}), \cdots, (\theta^{(S)}, \sigma^{2^{(S)}}) \}$` represent samples from the true posterior? As $S \rightarrow \infty$, this empirical distribution approaches the true posterior $p(\theta, \sigma^2 \mid y)$.

The Gibbs sampler is one example of a `Markov chain Monte Carlo` (MCMC) sampling method. The development of MCMC methods, combined with improved computing resources, has brought Bayesian statistics back into popularity. This will become clearer as we discuss various other MCMC methods in the coming lectures.

These results also suggest that one can approximate posterior summaries of interest using these samples, as long as $S$ is *large enough*. For example, the posterior mean for $\theta$ is

<div>
$$
    E(\theta \mid y) \approx \frac{1}{S} \sum_{s=1}^S \theta^{(s)}
$$
</div>

Similarly we can get a posterior mean for $\sigma^2$:

<div>
$$
    E(\sigma^2 \mid y) \approx \frac{1}{S} \sum_{s=1}^S \sigma^{2^{(s)}}
$$
</div>

The only thing that's different from before is the way we're generating those sequences.

## Back to the example

Going back to the cereals example, recall that our sampling model is

<div>
$$
    Y_1, \cdots, Y_n \mid \theta, \sigma^2 \overset{i.i.d.}{\sim} N(\theta, \sigma^2)
$$
</div>

And we specified independent priors:

<div>
$$
    \begin{gathered}
        \theta \sim N(\mu_0, \tau_0^2) \\
        \sigma^2 \sim IG(a, b)
    \end{gathered}
$$
</div>

### Prior hyperparameters selection

The hyperparameters that we have to choose are $\mu_0$, $\tau_0^2$, $a$ and $b$. Suppose we have little prior information about the parameters, except a guess that the average calories per serving is maybe near 200. This might suggest the following hyperparameter values:

- $\mu_0=200$, $\tau_0^2 = 65^2$ to ensure virtually zero probability of $\theta < 0$.
- $a=b=0.01$ as a common choice for the prior on $\sigma^2$ that has little influence.

### Full conditional distributions

{{<hl>}}The full conditional distribution is always proportional to the posterior distribution{{</hl>}}. The full conditional distribution for $\theta$ is the same as the posterior we got in the unknown $\theta$, fixed $\sigma^2$ case:

<div>
$$
    \theta \mid \sigma^2, y \sim N(\mu_n, \tau_n^2)
$$
</div>

where:

<div>
$$
    \mu_n = \cfrac{\cfrac{\mu_0}{\tau_0^2} + \cfrac{n\bar{y}}{\sigma^2}}{\cfrac{1}{\tau_0^2} + \cfrac{n}{\sigma^2}}, \quad \tau_n^2 = \cfrac{1}{\cfrac{1}{\tau_0^2} + \cfrac{n}{\sigma^2}}
$$
</div>

Plugging in $n=77$, $\bar{y} = 106.88$ and prior hyperparameters, we can find both $\mu_n$ and $\tau_n^2$ to be functions of $\sigma^2$.

Similarly, the full conditional distribution for $\sigma^2$ is similar to the posterior we got in the fixed $\theta$, unknown $\sigma^2$ case:

<div>
$$
    \sigma^2 \mid \theta, y \sim IG(a_n, b_n)
$$
</div>

where:

<div>
$$
    a_n = a + \frac{n}{2}, \quad b_n = b + \frac{1}{2}\sum_{i=1}^n (y_i - \theta)^2 = b + \frac{1}{2}\left[ (n-1)s^2 + n(\bar{y} - \theta)^2 \right]
$$
</div>

Plugging in the known values, both $a_n$ and $b_n$ are functions of $\theta$. Now for the Gibbs sampler, we can initialize the parameters[^init-choice] with

[^init-choice]: These are good choices as we're using the sample mean and variance to initialize the true mean and variance values. In practice these don't matter too much, especially if we draw a large number of samples. 

<div>
$$
    \theta^{(1)} = \bar{y}, \quad \sigma^{2^{(1)}} = s^2
$$
</div>

Then in each scan/sweep to update the parameters, we:

- Draw $\theta^{(s)} \sim N(\mu_n, \tau_n^2)$ where $\mu_n$ and $\tau_n^2$ are computed setting $\sigma^2 = \sigma^{2^{(s-1)}}$.
- Draw $\sigma^{2^{(s)}} \sim IG(a_n, b_n)$ where $a_n$ and $b_n$ are computed using $\theta = \theta^{(s)}$.
