---
title: "Metropolis-Hastings Algorithms"
date: 2021-03-22T10:20:00-05:00
summary: "This lecture discusses the Metropolis and Metropolis-Hastings algorithms, two more tools for sampling from the posterior distribution when we do not have it in closed form. These are used when we are unable to obtain full conditional distributions. MCMC for the win!" # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-metropolis-hastings-algorithms"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 80 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "San Francisco 2077." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "5CN_4tEHDw4" # Unsplash ID of the picture
---

This lecture discusses the Metropolis and Metropolis-Hastings algorithms, two more tools for sampling from the posterior distribution when we do not have it in closed form. However, unlike the Gibbs sampler, we use these methods when we are not able to obtain full conditional distributions. Discussion of intuition behind the algorithms, as well as choices of various settings to assist in algorithm convergence is included.

## Quick recap

In the previous section we had a general problem of not being able to write down the posterior distribution `$p(\theta \mid y)$` in closed form. With the Gibbs sampler, we can approximately generate samples from the posterior if the following items are true:

-   `$\theta = (\phi_1, \cdots, \phi_p)$` is multi-dimensional, i.e. `$p > 1$`.
-   _All_ full conditional distributions `$p(\phi_i \mid \phi_{-i}, y)$` must be available in closed-form, where `$\phi_{-i}$` means all parameters expect `$\phi_i$`.

However, not every posterior we get is going to satisfy these. What if `$\theta$` is a single parameter, or in more generality, not all full conditional distributions are available in closed-form?

## Motivating example

Suppose we are examining a manufacturing process for the production of extension cords, which are meant to be 10 feet long. To do this, `$n=28$` cords are randomly selected from the production line, and the _error in lengths_ compared to the truth are measured. For example, if a cord is 9.9 feet long, then the error is -0.1 feet. 

Suppose we also know from a prior test that:

1.  The true standard deviation in errors is 0.05.
2.  The process is very likely to be unbiased, i.e. the process isn't known to produce extension cords which are systematically longer or shorter than expected.

The data is given as follows:


```r
library(tidyverse)
library(ggpubr)

error <- c(
  0.059210131, -0.066246887, 0.038790819, 0.06652549,
  -0.003005251, -0.031189845, 0.016891891  , 0.042535987,
  -0.046366324, -0.023882362, -0.018044908 , 0.103928565, 
  0.035109448, -0.013638036, -0.045720695, 0.082342375,
  0.076147619, 0.06187386, 0.027775895  , 0.016442953,
  0.014963665, 0.014737389, 0.003129473 , 0.050911582,
  0.052515172, -0.044619958, -0.000357528, 0.134953716
)
```

The question is should the process continue as normal, or is there strong evidence to suggest it has become biased? The sample mean and standard deviation are `$\bar{y} = 0.0216$` ft and `$s = 0.0489$` ft, respectively.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/err-density-1.svg" alt="Kernel density estimate of the errors in length." width="480" />
<p class="caption">Figure 1: Kernel density estimate of the errors in length.</p>
</div>

## Bayesian ingredients

Let `$Y_i$` be the error in length (from 10 feet) of sampled extension cord `$i$`, for `$i = 1, \cdots, 28$`. Let `$\theta$` be the true mean error, and `$\sigma^2$` be the true variance in the errors.

### Sampling model

From the prior test, it is assumed that `$\sigma^2 = 0.05^2$` is known. Our sampling model is:

<div>
$$
    Y_1, \cdots, Y_n \mid \theta \overset{i.i.d.}{\sim} \mathcal{N}(\theta, \sigma^2 = 0.05^2)
$$
</div>

The joint density is:

<div>
$$
    p(y \mid \theta) = \prod_{i=1}^n p(y_i \mid \theta) = (2\pi \times 0.05^2)^\frac{n}{2} \exp \left\{ -\frac{1}{2(0.05)^2} \sum_{i=1}^n (y_i - \theta)^2 \right\}
$$
</div>

### Prior model

We also have very strong evidence from the prior test that this manufacturing process is unbiased. For conjugacy purposes, we could try `$\theta \sim \mathcal{N}(0, \tau^2)$` for a very small `$\tau$`, e.g. `$\tau = 0.05$`. We center the prior at zero because our prior knowledge (unbiased) suggests `$\theta = 0$`. Note that under this prior, slight biases (close to `$\theta = 0$`) have a density nearly as high as `$p(0)$`[^prior-density]. This indicates that the prior might not put enough emphasis on `$\theta = 0$`.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/prior-density-1.svg" alt="Slight biases in the normal density." width="480" />
<p class="caption">Figure 2: Slight biases in the normal density.</p>
</div>

[^prior-density]: R code for slight biases in the normal density:

    
    ```r
    dnorm_label <- tibble(
      x = c(0, 0.025),
      y = map_dbl(x, ~ dnorm(.x, sd = 0.05)),
      label = fct_inorder(str_glue("p({x})"))
    )
    
    ggplot() +
      geom_function(fun = ~ dnorm(.x, sd = 0.05)) +
      geom_point(aes(x, y, color = label), data = dnorm_label, shape = 4, size = 2) +
      geom_segment(
        aes(-0.2, y, xend = x, yend = y, color = label),
        data = dnorm_label, linetype = "dashed"
      ) +
      geom_segment(
        aes(x, 0, xend = x, yend = y, color = label),
        data = dnorm_label, linetype = "dashed"
      ) +
      xlim(-0.2, 0.2) +
      labs(x = expression(theta), y = "Density") +
      theme_pubr() +
      ggsci::scale_color_nejm()
    ```

We might prefer a prior on `$\theta$` which captures our very strong belief that the true bias is zero. One possibility as an alternative is the `Laplace distribution`:

$$
    \theta \sim \text{Laplace}(\mu = 0, b = 0.01),
$$

which has the probability density function

<div>
$$
    p(\theta) = \frac{1}{2b} \exp \left\{ -\frac{|\theta - \mu|}{b} \right\}, \quad -\infty < \theta < \infty
$$
</div>

The Laplace distribution has mean `$E(\theta) = \mu$` and variance `$Var(\theta) = 2b^2$`. This is also known as the `double-exponential distribution`, because its PDF is shaped like two exponential distributions glued together. The Laplace distribution is sometimes used as a prior for regression coefficients when we have a large number of covariates compared to observations.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/normal-vs-laplace-1.svg" alt="Comparison of the Laplace prior to the normal prior." width="480" />
<p class="caption">Figure 3: Comparison of the Laplace prior to the normal prior.</p>
</div>

Comparing the Laplace prior to the normal prior[^laplace-vs-normal], the Laplace prior has a more prominent spike at `$\theta = 0$`. The Laplace distribution has heavier tails than the normal distribution as it is expressed in terms of the absolute difference from the mean rather than the squared difference.

[^laplace-vs-normal]: R code for Laplace vs. normal:

    
    ```r
    ggplot() +
      geom_function(fun = ~ dnorm(.x, sd = 0.05), color = "#BC3C29") +
      geom_function(fun = ~ exp(-abs(.x) / 0.01) / 0.02, color = "#0072B5") +
      geom_text(
        aes(x = 0.08, y = 8, label = "N(0, 0.05^2)"), parse = T, color = "#BC3C29"
      ) +
      geom_text(
        aes(x = 0.08, y = 40, label = "Laplace(0, 0.01)"), color = "#0072B5"
      ) +
      xlim(-0.2, 0.2) +
      labs(x = expression(theta), y = "Density") +
      theme_pubr()
    ```

### Computing the posterior

By Bayes' Theorem, the posterior is given by:

<div>
$$
    \begin{aligned}
        p(\theta \mid y) &\propto p(y \mid \theta) p(\theta) \\
        &\propto \exp\left\{ -\frac{1}{2(0.05)^2} \sum_{i=1}^n (y_i - \theta)^2 \right\} \exp\left\{ -\frac{|\theta|}{0.01} \right\} \\
        &\propto \exp\left\{ -\frac{1}{2(0.05)^2} \sum_{i=1}^n (y_i - \theta)^2 - \frac{|\theta|}{0.01} \right\}
    \end{aligned} 
$$
</div>

This is not recognizable as a closed-form distribution. Our next step might be to think about Gibbs sampling, but when we have only one parameter, the full conditional distribution for `$\theta$` is the same as the posterior for `$\theta$`, which we do not have in closed-form.

How do we compute posterior summaries of interest if we only know the posterior up to a normalizing constant? By Bayes' Theorem, we have:

<div>
$$
    p(\theta \mid y) = \frac{p(y \mid \theta) p(\theta)}{p(y)} = \frac{u(\theta \mid y)}{p(y)}
$$
</div>

where `$u(\theta \mid y)$` is the `un-normalized posterior density` (what we have) and `$p(y)$` is the normalizing constant (something we can't compute directly). Since the shape of the posterior distribution only depends on `$u(\theta \mid y)$`, can we generate samples from `$p(\theta \mid y)$` only using `$u(\theta \mid y)$`?

## Metropolis algorithm

Our intuition is that we want to generate approximate samples `$\{ \theta^{(1)}, \cdots, \theta^{(S)} \}$` from the posterior via Monte Carlo sampling for some large `$S$`. Just like the Gibbs sampler, we will do this in a **dependent** way where `$\theta^{(s)}$` depends on its previous value `$\theta^{(s-1)}$`. This requires an initial value of the parameter `$\theta^{(1)}$`.

At step `$s$`, the Gibbs sampler selects the next value `$\theta^{(s+1)}$` from the full conditional distribution of `$\theta$`. However, we don't have the full conditional here. Suppose we only know our current state `$\theta^{(s)}$`, and say we have a "guess" for `$\theta^{(s+1)}$`, which we will call `$\theta^*$`. We then have two ways of setting the next value:

<div>
$$
    \theta^{(s+1)} = \begin{cases}
        \theta^*, &\text{guess}, \\
        \theta^{(s)}, &\text{current state}
    \end{cases}
$$
</div>

To evaluate which one to select, there's two scenarios to consider. If `$\theta^*$` has higher posterior density than `$\theta^{(s)}$`,

<div>
$$
    p(\theta^* \mid y) > p(\theta^{(s)} \mid y) \implies \frac{p(\theta^* \mid y)}{p(\theta^{(s)} \mid y)} > 1,
$$
</div>

then we should set `$\theta^{(s+1)} = \theta^*$` because the guess is better than the current value. In the other case where `$\theta^*$` has lower posterior density than `$\theta^{(s)}$`:

<div>
$$
    p(\theta^* \mid y) < p(\theta^{(s)} \mid y) \implies \frac{p(\theta^* \mid y)}{p(\theta^{(s)} \mid y)} < 1,
$$
</div>

If we simply set the next value to be `$\theta^{(s)}$`. then if we keep proposing guesses that are bad, we never evolve from the current state. Furthermore, if our current state is the peak of the posterior, then we never explore the posterior distribution because all samples would be at the peak.

In practice, if `$\theta^*$` is close to `$\theta^{(s)}$` in terms of the associated densities, then there's a chance we might want to move to `$\theta^*$`. However, if `$\theta^*$` has a much smaller density, then we probably don't want to move to this worse guess.

In conclusion, `$p(\theta^* \mid y) > p(\theta^{(s)} \mid y)$` means if we've deemed `$\theta^{(s)}$` to be a plausible sample from our posterior, `$\theta^*$` is more plausible and we should **always** move to this value. `$p(\theta^* \mid y) < p(\theta^{(s)} \mid y)$` means `$\theta^{(s)}$` is less plausible than our current state, but we should still move to it at a fraction of the time depending on the ratio

<div>
$$
    \frac{p(\theta^* \mid y)}{p(\theta^{(s)} \mid y)}
$$
</div>

This properly reflects that `$\theta^{(s)}$` might be an appropriate value to sample for the posterior distribution, but it's just less plausible than our current state.

### Acceptance ratio

With all this, see that the posterior ratio of interest can be evaluated without knowing `$p(y)$`:

<div>
$$
    \frac{p(\theta^* \mid y)}{p(\theta^{(s)} \mid y)} = \cfrac{\cfrac{u(\theta^* \mid y)}{p(y)}}{\cfrac{u(\theta^{(s)} \mid y)}{p(y)}} = \frac{u(\theta^* \mid y)}{u(\theta^{(s)} \mid y)}
$$
</div>

Thus we can consolidate all our intuition about the value of the new sample `$\theta^{(s+1)}$` given the current state `$\theta^{(s)}$` into the following `acceptance ratio` involving the un-normalized posterior density:

<div>
$$
    r = \frac{u(\theta^* \mid y)}{u(\theta^{(s)} \mid y)}
$$
</div>

If `$r>1$`, then we always set `$\theta^{(s+1)} = \theta^*$`. If `$r<1$`, then we set `$\theta^{(s+1)} = \theta^*$` with probability `$r$`. This is equivalent to:

<div>
$$
    \theta^{(s+1)} = \begin{cases}
        \theta^*, & \text{with probability } \min\{1, r\} \\
        \theta^{(s)}, & \text{with probability } 1 - \min\{1, r\}
    \end{cases}
$$
</div>

### Proposal distribution

So how do we propose a guess for the new value given the current state? The answer is we generate `$\theta^*$` from a `proposal distribution`, with PDF denoted `$q(\theta \mid \theta^{(s)})$`.

The proposal distribution is usually chosen to depend on `$\theta^{(s)}$` so that we propose values close to the current state. This distribution is **not necessarily related to the original model we specify**! Two typical simple distributions are:

-   `$\theta^* \mid \theta^{(s)} \sim \mathcal{U}(\theta^{(s)} - c, \theta^{(s)} + c)$` for some pre-selected value `$c$`, and
-   `$\theta^* \mid \theta^{(s)} \sim \mathcal{N}(\theta^{(s)}, c^2)$` for some pre-selected value `$c$`.

If `$q(\theta^* \mid \theta^{(s)}) = q(\theta^{(s)} \mid \theta^*)$`, then we call this a `symmetric proposal` since it is equally likely to move from `$\theta^{(s)}$` to `$\theta^*$` as it is to go in the opposite direction.

### Metropolis algorithm

Suppose we have only the un-normalized posterior `$u(\theta \mid y)$`, i.e.

<div>
$$
    p(\theta \mid y) \propto u(\theta \mid y) = p(y \mid \theta) p(\theta)
$$
</div>

Then, we use the following steps to generate approximate samples from the posterior:

1. Set initial parameter value `$\theta^{(1)}$`.
2. Suppose `$\theta^{(s)}$` is the current state. Draw `$\theta^*$` from the symmetric proposal distribution `$q(\theta \mid \theta^{(s)})$`.
3. Compute the acceptance ratio:

 <div>
 $$
     r = \frac{u(\theta^* \mid y)}{u(\theta^{(s)} \mid y)}
 $$
 </div>

4. Set `$\theta^{(s+1)} = \theta^*$` with probability `$\min\{1, r\}$`; else, set `$\theta^{(s+1)} = \theta^{(s)}$`.
5. Repeat steps 2-4 until `$S$ `samples have been obtained.

### Metropolis-Hastings algorithm

The only difference between the `Metropolis-Hastings` algorithm and the Metropolis algorithm is that the proposal distribution is not necessarily symmetric in Metropolis-Hastings. Without the symmetry, we have different probabilities moving from one direction to the other. As a consequence, we could over-sample a part of the posterior. To adjust for this, the acceptance ratio becomes:

<div>
$$
    r = \frac{u(\theta^* \mid y)}{u(\theta^{(s)} \mid y)} \cdot \frac{q(\theta^{(s)} \mid \theta^*)}{q(\theta^* \mid \theta^{(s)})}
$$
</div>

## Back to the example

Using the Metropolis algorithm with initial value `$\theta^{(1)} = 0$` and proposal distribution:

<div>
$$
    \theta^* \mid \theta^{(s)} \sim \mathcal{N}(\theta^{(s)}, c^2), \quad c = 0.05,
$$
</div>

we can obtain the posterior from `$S=10000$` samples using the R code below. The explanation of the parameters can be found in the next section.


```r
library(patchwork)

n <- length(error)  # sample size
sigma <- 0.05   # param for sampling model
b <- 0.01       # hyperparam for prior

# Parameters for running Metropolis-Hastings algorithm
burnin_steps <- 0         # num of steps to "burn-in" each chain
n_chains <- 3             # num of chains to run
num_saved_steps <- 10000  # total num of steps to save
thin_steps <- 1           # num of steps to "thin" (1=keep every step)
n_iter <- ceiling((num_saved_steps * thin_steps) / n_chains)  # steps / chain
c <- 0.05                 # sd of proposal distribution (assuming normal)

# Initialize vector of sampled unknown parameter values, and select starting
# value for each chain. R is column-major so save each chain in a column
theta <- matrix(0, nrow = burnin_steps + n_iter, ncol = n_chains)
theta[1, ] <- rep(0, n_chains)   # initializes each chain to 0 currently

# Generate samples via Metropolis algorithm
accept <- rep(0, n_chains)      # num of accepted proposals

for (s in seq(2, burnin_steps + n_iter)) {
  for (i in seq(n_chains)) {
    # Propose a new value for theta (symmetric normal proposal distribution)
    theta_star <- rnorm(1, theta[s-1, i], c)

    # Compute acceptance ratio using logarithms - more stable numerically
    r <- exp(
      -(1 / (2*sigma^2)) * sum((error - theta_star)^2) -
      abs(theta_star) / b +
      (1 / (2*sigma^2)) * sum((error - theta[s-1, i])^2) +
      abs(theta[s-1, i]) / b
    )

    # Accept theta_star with probability min{1,r}, otherwise reject
    # Equivalent to flipping a coin with probability min{1,r}
    if (rbinom(1, 1, min(1, r)) == 1) {
      theta[s, i] <- theta_star
      if (s > burnin_steps) {
        accept[i] <- accept[i] + 1     # add one to "acceptance counter"
      }
    } else {
      theta[s, i] <- theta[s-1, i]
    }
  }
}

# Form posterior samples by removing burn-in period and thinning
keep_vals <- seq(burnin_steps + 1,burnin_steps + n_iter, by = thin_steps)
theta_post <- theta[keep_vals, ]

# Show trace plots of all chains simultaneously with posterior density
res <- theta_post %>%
  as.data.frame() %>%
  mutate(Iteration = keep_vals) %>%
  pivot_longer(-Iteration, names_to = "Chain", values_to = "theta") %>%
  mutate(Chain = str_remove(Chain, "^V"))

p1 <- ggline(res, x = "Iteration", y = "theta", color = "Chain",
             numeric.x.axis = T, plot_type = "l", palette = "nejm",
             ylab = expression(theta))

p2 <- ggdensity(res, x = "theta", add = "mean",
                xlab = expression(theta), ylab = "Density")
p1 | p2
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/metropolis-algo-1.svg" alt="Trace plots and posterior density of 10000 samples obtained using the Metropolis algorithm." width="960" />
<p class="caption">Figure 4: Trace plots and posterior density of 10000 samples obtained using the Metropolis algorithm.</p>
</div>

The final posterior density doesn't look exactly like a normal distribution because of the Laplace prior. We can use the samples in the posterior to compute posterior summaries:

```r
mean(theta_post)
# Posterior mean: 0.01373471

quantile(theta_post, c(0.025, 0.975))
# 95% credible interval:
#         2.5%        97.5%
# -0.001995037  0.032360048
```

The posterior mean is close to the sample mean. The 95% CI includes zero, so maybe this is not strong enough evidence to suggest that the process is biased.

## General considerations for MCMC algorithms

The Gibbs sampler, Metropolis, and Metropolis-Hastings algorithms are _guaranteed_ to sample from the true posterior distribution as `$S \rightarrow \infty$`. In practice, we only run the algorithm for a finite `$S$`. So how do we ensure our final collection of samples constitutes draws from the true posterior?

### Burn-in samples

Note that initial samples may not be reflective of the true posterior! It's typically recommended to remove a pre-specified number of initial samples, known as the `burn-in samples`, with the belief that the remaining samples are approximately settled around the true posterior. If the initializing parameters are far from the truth, then it may require longer to burn-in.

We can diagnose convergence via trace plots. The trace plot shown in the example above is what we want to see -- the random oscillations with no consistent up or down drifts. This indicates the Metropolis algorithm has converged/targeted the true posterior correctly.

If we initialized by a different value, e.g. `$\theta^{(1)} = 100$`, the trace plot looks bad. We have a considerable amount of drift downward and a tail in the posterior of `$\theta$` that goes to 100.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/metropolis-bat-init-1.svg" alt="Trace plot and posterior density with bad initialization." width="960" />
<p class="caption">Figure 5: Trace plot and posterior density with bad initialization.</p>
</div>

This happened because the initial value is far from values of `$\theta$` which we believe are likely to have posterior density. Our proposal distribution compounds that -- guesses with variance of `$0.05^2$` means we're moving very slowly away from 100. It takes the Metropolis algorithm a lot of time to settle down and sample from the true posterior, although it eventually converged.

If we separate the trace plots for the first 5500 iterations and the remaining 4500 iterations separately, the latter would look like random fluctuations. If we had specified the first 5500 iterations as burn-in, the algorithm would also appear to has converged.

There is no commonly accepted amount of burn-in to specify, and we usually have to determine the value based off the trace plots. It's okay to re-run MCMC multiple times, because we're not changing the model itself -- we're only changing how we sample from the model.

### Proposal distribution selection

The proposal distribution used for the manufacturing example is:

<div>
$$
    \theta^* \mid \theta^{(s)} \sim \mathcal{N}(\theta^{(s)}, c^2)
$$
</div>

Setting `$c$` too small means the proposed values will be very close to the current state. Most proposed values will be accepted, but the trade-off is the exploration of the posterior is very slow.

Setting `$c$` too large means proposed values will be very far from the current state. If the current state is in a region of plausible values, the proposal may jump to a value with low posterior density. Most proposed values will not be accepted, and that again results in very slow exploration of the posterior.

To address this problem, we may:

-   `Thinning` the sample. For example, we only take every 5th sample as a member of the approximate posterior. This helps when `$c$` is too small, because the guesses depend a lot on the current state. This can help ensure the posterior samples are "independently drawn".
-   Track the `percentage of accepted proposals`. In general, acceptance rate will be very high when `$c$` is too small, and low when `$c$` is too large. Some have shown that an acceptance rate close to **23.4%** is optimal for convergence.

    `Adaptive MCMC algorithms` let the value `$c$` update to target this optimal acceptance rate. One example of such software is called JAGS (Just Another Gibbs Sampler).

### Diagnose convergence

When we are trying to diagnose the convergence of MCMC samplers, one thing that's commonly done in practice is to run multiple chains independently, initializing each one at a different `$\theta^{(1)}$`. If the samples have truly converged to the posterior distribution, then the trace plots for all chains should coincide with each other[^gelman-rubin].

In general, we want to run the MCMC algorithm for as long as possible, and also run several chains with different initial values. We may never know if MCMC has converged to the correct posterior, but we will definitely know if something is wrong!

[^gelman-rubin]: This can be put into a hypothesis test framework -- the Gelman-Rubin diagnostic statistic compares between-chain to within-chain variance. Not many people seem to use it though.

### Metropolis vs. Metropolis-Hastings

As we stated earlier, Metropolis requires a symmetric proposal distribution, whereas Metropolis-Hastings can be used for asymmetric proposal distributions. Metropolis is a special case of Metropolis-Hastings.

So why would one bother with asymmetric proposal distributions? One case is when the parameter space is bounded, e.g. `$\theta$` must be between 0 and 1. If our current state is close to 1, then we don't want to propose a value that's greater than 1. In this case, we might favor moving towards smaller values of `$\theta$` to get away from the boundary of the parameter space. Introducing asymmetry can also improve how quickly the chain correctly samples from the true posterior (improved convergence rate).

Some modern samplers which use asymmetric proposal distributions are Metropolis-adjusted Langevin algorighm (MALA), Hamiltonian Monte Carlo (HMC), and reversible-jump MCMC. These often come from physics applications such as energy function optimization. They are often more computationally efficient, but are also more complicated in terms of how we propose things.
