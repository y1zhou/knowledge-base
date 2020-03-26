---
title: "Bias and Variance"
slug: "mathematical-statistics-bias-and-variance"
categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Estimation

summary: ""
date: 2020-01-25T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 110

menu:
  maths-stat:
    name: Bias and Variance
    parent: Estimation
    weight: 10
---



{{% alert note %}}

To get started, we're going to introduce some basic concepts.

{{% /alert %}}



## Basic concepts

**Def:** A `parameter` is a numerical characteristic of a population distribution, and is often *unknown*. A `Statistic` is a numerical summary of the sample, and it should depend on the sample only and does **not** involve unknown parameters.

If we have a model $Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} Bern(p)$, here $p$ is unknown and is a parameter. $Y_1$ is a statistic of the sample, and so is $Y_1 + \cdots + Y_n$. However, $\sum_{i=1}^n (Y_i - p)^2$ is not a statistic because the unknown parameter $p$ is involved.

More rigorously speaking, a quantity $T$ is called a `statistic` if it can be expressed as a function of the sample

$$
T = f(Y_1, \cdots, Y_n)
$$

where $f$ is known. {{<hl>}}The sample size $n$ is always treated as known.{{</hl>}} A statistic is a fixed rule to calculate a quantity based on the sample. The same sample values yields the same statistic values.

**Def:** `Point estimation` is a single value estimate of a parameter $\theta$ based on the sample. Often, an estimate of $\theta$ is called an `estimator` and is denoted $\hat{\theta}$. An estimator is a statistic as it doesn't involve anything unknown. Because of this, an estimator is also a function of the sample.

**Example:** suppose we have an i.i.d. sample $Y_1, \cdots, Y_n$ following the same distribution with mean $\mu$ and variance $\sigma^2$. Are the following estimators?

| Estimator                                          | Known $\mu$ | Unknown $\mu$ |
| -------------------------------------------------- | ----------- | ------------- |
| $T_1 = \frac{1}{n} \sum_{i=1}^n (Y_i - \bar{Y})^2$ | Yes         | Yes           |
| $T_2 = \frac{1}{n} \sum_{i=1}^n (Y_i - \mu)^2$     | Yes         | No            |

A commonly used type of estimator is called an `interval estimate`. In general, this is an interval constructed based on the sample, which hopefully contains the true parameter with certain quantified accuracy. The lower and upper limits $L$ and $U$ do **not** involve unknown parameters.

**Example:** confidence intervals (a frequentist method where we assume the unknown parameters are fixed) and credible intervals (a Bayesian method where the unknown parameters are random).



---



Some notations before we start:

- $\theta$: unknown parameter
  - Taking values within a set $\Theta$ called the `parameter space`, which is often a subset of the real line $\mathbb{R}$.
  - E.g. $p$ in $Bern(p)$ has a parameter space of $\Theta = [0, 1]$.
- $\hat{\theta}$: an estimator of $\theta$. It has to be a known function of sample $(Y_1, \cdots, Y_n)$.

Now, a very basic question is how do we evaluate an estimator? Suppose we have an observed (realized) sample $(y_1, \cdots, y_n)$. We can plug these values into two estimators

$$
\begin{gather*}
	\hat{\theta}_1 = f_1(Y_1, \cdots, Y_n) \\\\
	\hat{\theta}_2 = f_2(Y_1, \cdots, Y_n)
\end{gather*}
$$

and get two point estimates. Let's say we know the true value of $\theta$ and the value of $\hat{\theta}_1$ was closer to it. Is $\hat{\theta}_1$ better than $\hat{\theta}_2$? The answer's no because we cannot evaluate an estimator based on a single realized estimate. However, $\hat{\theta}$ is a function of *random* sample $(Y_1, \cdots, Y_n)$. Some realizations may make $\hat{\theta}$ closer to $\theta$, while others may not.

To address the randomness, an idea is to repeat realizations of sample $(Y_1, \cdots, Y_n)$ and plot a histogram of the realized estimator. Comparing the histograms of different estimators may give us some insight on which is a better one.

## Bias

**Def:** The `bias` of estimator $\hat{\theta}$ of $\theta$ is defined as

$$
Bias(\hat{\theta}; \theta) = E[\hat{\theta}] - \theta,
$$

where the expectation $E$ is taken by assuming that $\theta$ is the true parameter ($E_{\theta}$).

**Remark:** the bias is a function of the parameter $\theta \in \Theta$. We'll come back to this when introducing the concept of unbiasedness.

### Uniform distribution example

$(Y_1, \cdots, Y_n) \overset{i.i.d.}{\sim} Uniform(0, \theta)$. We want to estimate $\theta$. The proposed estimator is $\bar{Y}_n$. The bias of our estimator is

$$
\begin{aligned}
	Bias(\bar{Y}_n; \theta) &= E[\bar{Y}_n] - \theta \\\\
	&= E\left[\frac{1}{n} (Y_1 + \cdots + Y_n)\right] - \theta \\\\
	&= \frac{1}{n} E[Y_1 + \cdots + Y_n] - \theta \\\\
	&= \frac{1}{n} \sum\_{i=1}^n E[Y_i] - \theta \\\\
	&= \frac{1}{n} \cdot n \cdot \frac{\theta}{2} - \theta \\\\
	&= -\frac{\theta}{2}
\end{aligned}
$$

Note that $E[\bar{Y}_n] = E[Y_i]$. This holds true when the samples are i.i.d. Note also that the bias is a function of $\theta$. Our parameter space is $\Theta = (0, \infty)$.

**Def:** An estimator $\hat{\theta}$ of $\theta$ is said to be `unbiased` if $Bias(\hat{\theta}; \theta) = 0$, i.e. $E[\hat{\theta}] = \theta$ for any $\theta \in \Theta$. The "for any" part is important because $\theta$ is unknown.

### Unbiasedness example

Suppose $(Y_1, \cdots, Y_n)$ is i.i.d. and $E[Y_i] = \theta \in \Theta = \mathbb{R}$. $\bar{Y}_n$ is an estimator of $\theta$. We always have $E[\bar{Y}_n] = \theta$, so the bias of this estimator is $0$ for all $\theta \in \mathbb{R}$. According to the definition, the estimator $\bar{Y}_n$ is unbiased for $\theta$.

{{<hl>}}Unbiasedness is a desirable property. It says, averagely speaking, the estimator captures the true parameter no matter where the latter is located.{{</hl>}}

### Biased example

$(Y_1, \cdots, Y_n) \overset{i.i.d.}{\sim} Unif(0, \theta)$. The proposed estimator $\hat{\theta}_1 = 1$. Suppose the true $\theta = 1$. Then the bias would be

$$
Bias(\hat{\theta}_1; \theta) = E[\hat{\theta}_1] - \theta = 1 - 1 = 0.
$$

This is **wrong** because $\theta$ is unknown and we can't just assume its value. What we know is that $\Theta = (0, \infty)$. Suppose $\theta \in \Theta$,

$$
Bias(\hat{\theta}_1; \theta) = E[\hat{\theta}_1] - \theta = 1 - \theta.
$$

So the conclusion is $\hat{\theta}_1$ is biased. Here we emphasize the importance of the "for all" language in the definition of unbiasedness.

Another proposed estimator is $\hat{\theta}_2 = 2\bar{Y}_n$. Since we know $E[\bar{Y}_n] = \frac{\theta}{2}$,

$$
E[\hat{\theta}_2] = 2E[\bar{Y}_n] = \theta.
$$

Thus, $Bias(\hat{\theta}_2; \theta) \equiv 0$ and $\hat{\theta}_2$ is unbiased.

**Proposition:** suppose $\hat{\theta}$ is an unbiased estimator for $\theta$, and $\hat{\eta}$ is an unbiased estimator for $\eta$.

1. Let $a$ be a known non-random scalar. $a\hat{\theta}$ would be unbiased for $a\theta$.
2. $\hat{\theta} + a$ would be unbiased for $\theta + a$.
3. $a\hat{\theta} + b\hat{\eta}$ would be unbiased for $a\theta + b\eta$.

In conclusion, {{<hl>}}unbiasedness is preserved under linear transformations. However, it's often **not** preserved under nonlinear transformations.{{</hl>}}

### Nonlinear transformation example

Suppose $\hat{\theta}$ is an unbiased estimator for $\theta$. Is $\hat{\theta}^2$ an unbiased estimator of $\theta^2$?

$$
E\left[\hat\theta^2\right] = Var(\hat\theta) + \left( E\left[\hat\theta\right] \right)^2 =  Var(\hat\theta) + \theta^2
$$

Unless $ Var(\hat\theta) = 0$, we always have a bias. $ Var(\hat\theta) = 0$ implies $\hat\theta$ is a constant, which is very rare.

**Proposition:** $Y_1, \cdots, Y_n$ ($n \geq 2$) are i.i.d. random samples with population variance $\sigma^2$. The estimator
$$
\hat{\sigma}_n^2 = \frac{1}{n-1}\sum\_{i=1}^n (Y_i - \bar{Y}_n)^2
$$

is an unbiased estimator of $\sigma^2$. The proof is given as follows.

$$
\begin{aligned}
	S &= \sum\_{i=1}^n \left(Y_i - \bar{Y}\right)^2 \\\\
	&= \sum\_{i=1}^n \left(Y_i^2 - 2Y_i\bar{Y} + \bar{Y}^2\right) \\\\
	&= \sum\_{i=1}^n Y_i^2 - 2\bar{Y}\underbrace{\sum_{i=1}^n Y_i}_{n\bar{Y}} + n\bar{Y}^2 \\\\
	&= \sum\_{i=1}^n Y_i^2 - n\bar{Y}^2 = A - B
\end{aligned}
$$

Recall that $E\left[ Y_i^2\right] = \sigma^2 + \mu^2$ where $\mu = E[Y_i]$.

$$
\begin{aligned}
	E[A] &= \sum\_{i=1}^n E[Y_i^2] = n\sigma^2 + n\mu^2 \\\\
	Var\left(\bar{Y}\right) &= Var\left( \frac{Y_1 + \cdots + Y_n}{n} \right) \\\\
	&= \frac{1}{n^2}\left( Var(Y_1) + \cdots + Var(Y_n) \right) \\\\
	&= \frac{1}{n^2}n\sigma^2 = \frac{\sigma^2}{n} \\\\
	E\left[ \bar{Y}^2 \right] &= Var(\bar{Y}) + E[\bar{Y}]^2 = \frac{\sigma^2}{n} + \mu^2 \\\\
	\Rightarrow E[B] &= \sigma^2 + n\mu^2 \\\\
	E\left[ \frac{S}{n-1} \right] &= \frac{1}{n-1}(E[A] - E[B]) \\\\
	&= \sigma^2
\end{aligned}
$$

## Variance

The importance of the variance is best explained with an example. Suppose $Y_1, \cdots, Y_n$ are i.i.d. with mean $\mu$ and variance $\sigma^2$. Previously we've used $\bar{Y}$ to estimate $\mu$. Why didn't we use $Y_1$ if $E[Y_1] = \mu$? Well,

$$
Var(\bar{Y}) = \frac{\sigma^2}{n} \text{ vs. }Var(Y_1) = \sigma^2
$$

We can easily see that when $n > 2$, $Var(\bar{Y}) < Var(Y_1)$.

**Def:** The `variance` of estimator $\hat\theta$ of $\theta$ is

$$
Var(\hat\theta; \theta) = Var(\hat\theta),
$$

where the variance is computed by assuming that $\theta$ is the true parameter. It is also a function of $\theta$. The square root of the variance is denoted

$$
S.E.(\hat\theta) = \sqrt{Var(\hat\theta; \theta)}
$$



{{% alert note %}}

In this and later chapters, "standard error" and "standard deviation" are used interchangeably.

{{% /alert %}}

### Bernoulli example

$Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} Bern(p)$. We know that $E[Y_i] = p$ and $Var(Y_i) = p(1-p)$. Let $X_n = Y_1 + \cdots + Y_n$. The estimator

$$
\hat{p} = \bar{Y} = \frac{X_n}{n}
$$

is unbiased for $p$. Find the variance of $\hat{p}$.

$$
Var(\hat{p}) = Var(\bar{Y}) = \frac{Var(Y_1)}{n} = \frac{p(1-p)}{n}
$$

### Linear combination example

$\\{X_1, \cdots, X_{n_1}\\}$ are i.i.d. with $E[X_i] = \mu_1$ and $Var(X_i) = \sigma_1^2$. $\\{Y_1, \cdots, Y_{n_2}\\}$ are i.i.d. with $E[Y_i] = \mu_2$ and $Var(Y_i) = \sigma_2^2$. We assume that $\{X_i\}$ and $\{Y_i\}$ are independent of each other. We want to estimate $\theta = \mu_1 - \mu_2$.

Our proposed estimator is $\hat\theta = \bar{X} - \bar{Y}$. We know that $\bar{X}$ is unbiased for $\mu_1$ and $\bar{Y}$ is unbiased for $\mu_2$. Their linear combination would also be unbiased.

We're also interested in the variability of this estimator. We know that

$$
Var(A + B) = Var(A) + Var(B) + 2Cov(A, B),
$$

so the variance of $\hat\theta$ would be

$$
\begin{aligned}
	Var(\hat\theta) &= Var(\bar{X} - \bar{Y}) \\\\
	&= Var(\bar{X}) + Var(-\bar{Y}) + 2Cov(\bar{X}, \bar{Y}) \\\\
	&= \frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2} + 0 \\\\
	&= \frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}
\end{aligned}
$$

So the standard error of $\hat\theta$ is

$$
S.E.(\hat\theta) = \sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}}
$$

### Uniform distribution example

$Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} Unif(0, \theta)$. The two proposed estimators are $\hat{\theta}_1 = 1$ and $\hat{\theta}_2 = 2\bar{Y}$. Find the variance of the two estimators.

$$
\begin{gather*}
	Var(\hat\theta_1) = 0 \\\\
	Var(\hat\theta_2) = 4Var(\bar{Y}) = \frac{4Var(Y_1)}{n} = \frac{4\sigma^2}{12n} = \frac{\theta^2}{3n}
\end{gather*}
$$

We've shown earlier that $\hat\theta_1$ is biased, but here in terms of variability it's actually perfect. How do we determine which one is the better estimator?

## Mean squared error

**Def:** We need a measure of goodness for estimators combining both bias and variance. The `MSE` is defined as

$$
MSE(\hat\theta; \theta) = E\left[(\hat\theta - \theta)^2\right]
$$

where the expectation is taken by assuming that $\theta$ is the true parameter.

**Theorem:** $MSE(\hat\theta; \theta) = Bias(\hat\theta; \theta)^2 + Var(\hat\theta; \theta)$. This decomposition is derived as follows.

$$
\begin{aligned}
	MSE(\hat\theta) &= E\left[ \left((\hat\theta - E[\hat\theta]) + Bias(\hat\theta)\right)^2 \right] \\\\
    &= E\left[ (\hat\theta - E[\hat\theta])^2 \right] + 2E\left[ Bias(\hat\theta)(\hat\theta - E[\hat\theta]) \right] + E\left[ Bias(\hat\theta)^2 \right] \\\\
    &= Var(\hat\theta) + 2Bias(\hat\theta)\underbrace{E\left[\hat\theta - E[\hat\theta]\right]}_{0} + Bias(\hat\theta)^2 \\\\
    &= Var(\hat\theta) + Bias(\hat\theta)^2
\end{aligned}
$$

If $\hat\theta$ is unbiased,

$$
MSE(\hat\theta; \theta) = Var(\hat\theta; \theta)
$$

> **Notation:** for simplicity, we can write $Bias(\hat\theta)$, $Var(\hat\theta)$ and $MSE(\hat\theta)$ if it's clear from the context what's the assumed true parameter $\theta$.

A question here is why do we take the square? Why not $E\left[|\hat\theta - \theta|\right]$?

- The main reason is just mathematical convenience.
- Another reason is "[risk aversion](https://en.wikipedia.org/wiki/Risk_aversion)". The square penalizes large deviations more than the absolute value.

However, the other choices of "penalization functions" are **not** useless. These are further studied in robust statistics. Another note here is there's almost always a trade-off between the variance and the bias. We'll discuss this in detail later.

### Calculation example

Suppose $Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$. Three estimators are proposed:

1. $\hat\theta_1 = 1$,
2. $\hat\theta_2 = 2\bar{Y}_n$, and
3. $\hat\theta_3 = \max(Y_1, \cdots, Y_n)$.

Calculate the MSE of each estimator.

$$
\begin{gather*}
	MSE(\hat\theta_1) = Bias(\hat\theta_1)^2 + Var(\hat\theta_1) = (1 - \theta)^2 + 0 \\\\
	MSE(\hat\theta_2) = 0^2 + \frac{\theta^2}{3n} = \frac{\theta^2}{3n} \quad \cdots \rightarrow 0 \text{ as } n \rightarrow \infty
\end{gather*}
$$

In addition to $\theta$, the sample size also enters in the $MSE$ for $\hat\theta_2$. when $n$ is large, the $MSE$ becomes really small. Usually we fix $n$ at some value.

For $\hat\theta_3$, we're calculating the `1st order statistic`, which is often denoted $Y_{(1)}$. The first step is to compute the PDF of $\hat\theta_3$, which is also the hardest part of the problem. We know for a fact that

$$
\max(Y_1, \cdots, Y_n) \leq x \Leftrightarrow Y_1 \leq x, Y_2 \leq x, \cdots, Y_n \leq x
$$

and this is very useful in getting the CDF.
$$
\begin{aligned}
	P(\hat\theta_3 \leq x) &= P\Big(\max(Y_1, \cdots, Y_n) \leq x \Big), \quad 0 < x < \theta \\\\
	&= P(Y_1 \leq x, \cdots, Y_n \leq x) \\\\
	&= P(\{Y_1 \leq x\} \cap \cdots \cap \{Y_n \leq x\}) \\\\
	&= P(Y_1 \leq x) \cdot P(Y_2 \leq x) \cdots P(Y_n \leq x) \quad \cdots \text{independence} \\\\
	&= P(Y_1 \leq x)^n \\\\
	&= \left(\frac{x}{\theta}\right)^n
\end{aligned}
$$

By differentiation, the PDF is

$$
f(x) = \frac{d}{dx}\left(\frac{x}{\theta}\right)^n  =\frac{nx^{n-1}}{\theta^n}, \quad 0 < x < \theta
$$

Now we can find the bias and the variance.

$$
\begin{aligned}
	E[\hat\theta_3] &= \int_0^\theta xf(x)dx \\\\
	&= \int_0^\theta x \cdot \frac{nx^{n-1}}{\theta^n} dx \\\\
	&= \frac{n}{\theta^n} \cdot \frac{x^{n+1}}{n+1} \bigg|_0^\theta \\\\
	&= \frac{n}{n+1}\theta \\\\
	Bias(\hat\theta_3) &= -\frac{1}{n+1}\theta
\end{aligned}
$$

**Remark:** although $\hat\theta_3$ is biased, bias $\rightarrow 0$ as $n \rightarrow \infty$. We call this `asymtotically unbiased`.

To find $Var(\hat\theta_3)$, we need to find $E\left[ \hat\theta_3^2 \right]$ first.
$$
\begin{aligned}
	E\left[ \hat\theta_3^2 \right] &= \int_0^\theta x^2 \cdot \frac{nx^{n-1}}{\theta^n} dx \\\\
	&= \frac{n}{\theta^n} \int_0^\theta x^{n+1}dx \\\\
	&= \frac{n\theta^2}{n+2}
\end{aligned}
$$

Now we have


$$
\begin{aligned}
    Var(\hat\theta_3) &= E\left[ \hat\theta_3^2 \right] - E[\hat\theta_3]^2 \\\\
	  &= \frac{n\theta^2}{n+2} - \frac{n^2\theta^2}{(n+1)^2} \\\\
	&= \frac{n\theta^2}{(n+1)^2(n+2)}
\end{aligned}
$$


And now we can find the MSE:
$$
\begin{aligned}
	MSE(\hat\theta_3) &= Bias(\hat\theta_3)^2 + Var(\hat\theta_3) \\\\
	&= \frac{\theta^2}{(n+1)^2} + \frac{n\theta^2}{(n+1)^2(n+2)} \\\\
	&= \frac{2\theta^2}{(n+1)(n+2)}
\end{aligned}
$$

We can also find that

$$
\frac{MSE(\hat\theta_3)}{MSE(\hat\theta_2)} = \frac{6n}{(n+1)(n+2)} < 1 \text{ if } n > 2,
$$

which means $MSE(\hat\theta_3) < MSE(\hat\theta_2)$ when $n > 2$.

### Efficiency
**Def:** Given two estimators $\hat\theta_1$ and $\hat\theta_2$ of the same parameter based on the sample random sample, the `efficiency` of $\hat\theta_1$ relative to $\hat\theta_2$ is defined as
$$
eff(\hat\theta_1, \hat\theta_2) = \frac{MSE(\hat\theta_2)}{MSE(\hat\theta_1)}
$$

This works for both biased and unbiased estimators, although we usually look at the efficiency for unbiased estimators.
