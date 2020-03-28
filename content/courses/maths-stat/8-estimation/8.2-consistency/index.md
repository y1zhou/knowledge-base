---
title: "Consistency"
slug: "mathematical-statistics-consistency"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Estimation

summary: "Introducing consistency, a concept about the convergence of estimators. We start from the convergence of non-random number sequences to convergence in probability, then to consistency of estimators and its properties."
date: 2020-01-27T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 120

menu:
  maths-stat:
    name: Consistency
    parent: Estimation
    weight: 20
---



## Convergence

Consistency is about the convergence of estimators. Recall what convergence means for non-random numbers. Suppose $x_1, x_2, x_3, \cdots$ are non-random numbers. What is the meaning of $\lim\limits_{n \rightarrow \infty} x_n = x$?

For example, if $x_n = \frac{1}{n}$, $\lim\limits_{n \rightarrow \infty} x_n = 0$. If $x_n = (-1)^n$, $\lim\limits_{n \rightarrow \infty} x_n$ doesn't exist.

**Def:** A sequence $x_n \in \mathbb{R}$ is said to `converge` to $x$, denoted as $\lim\limits_{n \rightarrow \infty} x_n = x$ if for any fixed number $\epsilon > 0$, we have $|x_n - x| \leq \epsilon$ for all sufficiently large $n$. "For all sufficiently large $n$" means that there exists $N$ such that for all $n \geq N$.

### Example

Suppose $Y_i \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$. $\bar{Y}$ converges to $\mu$. We claim that $\bar{Y}_n \sim N(\mu, \frac{\sigma^2}{n})$. See [this theorem]({{< ref "/courses/maths-stat/5-functions-of-random-variables/index.md#theorem" >}}) for details.

No matter how large $n$ is, $\bar{Y}_n$ has a positive probability to exceed any fixed threshold (think about the bell-shaped curve). The good news is that

$$
\underbrace{MSE(\bar{Y}_n) = Var(\bar{Y}_n)}\_{\text{unbiasedness}} = \frac{\sigma^2}{n} \rightarrow 0 \text{ as } n \rightarrow \infty.
$$

## Convergence in probability

**Def:** A sequence of random variables $X_n$ is said to `converge in probability` to a constant $x$ if for *any* fixed $\epsilon > 0$,
$$
P\left( |X_n - x| \leq \epsilon \right) \rightarrow 1 \text{ as } n \rightarrow \infty.
$$

This is the same as

$$
P(x - \epsilon \leq X_n \leq x + \epsilon) \rightarrow 1,
$$

or

$$
P\left( |X_n - X| > \epsilon \right) \rightarrow 0 \text{ as } n \rightarrow \infty.
$$

The above (converge in probability to) can be denoted as

$$
X_n \xrightarrow{P} x
$$

{{% alert note %}}

The concept also applies to the case where the limit is also random. For us, the limit is always non-random in this course.

{{% /alert %}}



## Consistency

**Def:** An estimator $\hat\theta_n$ is said to be `consistent` if $\hat\theta_n \xrightarrow{P} \theta$ no matter which true $\theta \in \Theta$ is. Here $n$ denotes the sample size.

### Normal distribution example

Suppose $Y_i \overset{i.i.d.}{\sim}N(\mu, \sigma^2)$. We know that $\bar{Y} \sim N(\mu, \frac{\sigma^2}{n})$. Show the consistency of $\bar{Y}$.

$$
P(|\bar{Y}_n - \mu| \leq \epsilon) = P\left( \left| \frac{\bar{Y}_n - \mu}{\sigma / \sqrt{n}} \right| \leq \frac{\epsilon}{\sigma / \sqrt{n}} \right)
$$

The above is the `standardization` of a random variable. A fact here is that $Z = \frac{\bar{Y}_n - \mu}{\sigma / \sqrt{n}} \sim N(0, 1)$ because $E[\bar{Y}_n] = \mu$ and $s.e.(\bar{Y}_n) = \sigma / \sqrt{n}$.

$$
P(|\bar{Y}_n - \mu| \leq \epsilon) = P\left( |Z| \leq \frac{\epsilon}{\sigma/\sqrt{n}} \right) = P\left( |Z| \leq \frac{\epsilon}{\sigma}\sqrt{n} \right)
$$

which is the area under $\left( -\frac{\epsilon}{\sigma}\sqrt{n}, \frac{\epsilon}{\sigma}\sqrt{n} \right)$ in the PDF of $N(0, 1)$. As $n \rightarrow \infty$, the boundaries get pushed further outside and eventually we get the area to be $1$. 

### Uniform distribution example

$Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$. $\hat\theta_n = \max(Y_1, \cdots, Y_n)$. We want to show $\hat\theta_n$ is consistent.

We know that $\hat\theta_n \in [0, \theta]$. We can assume that $\epsilon \in (0, \theta)$ since $|\hat\theta_n - \theta| > \theta$ is impossible.

$$
\begin{aligned}
	P(|\hat\theta_n - \theta| > \epsilon) &= P(\theta - \hat\theta_n > \epsilon) \\\\
	&= P(\hat\theta_n < \theta - \epsilon) \\\\
	&= P(Y_1 < \theta - \epsilon, Y_2 < \theta - \epsilon, \cdots, Y_n < \theta - \epsilon) \\\\
	&= P(Y_1 < \theta - \epsilon)^n \\\\
	F_{Y_1}(y) &= \frac{y}{\theta}, 0 \leq y \leq 1 \\\\
	P(|\hat\theta_n - \theta| > \epsilon) &= \left( \frac{\theta - \epsilon}{\theta} \right)^n \rightarrow 0 \text{ as } n \rightarrow \infty.
\end{aligned}
$$

###  Theorem
If $MSE(\hat\theta_n; \theta) \rightarrow 0$ as $n \rightarrow \infty$ $\forall \theta \in \Theta$, then $\hat\theta_n$ is consistent.

**Lemma (Markov inequality):** If random variable $X \geq 0$, then for any constant $k > 0$, we have
$$
P(X > k) \leq \frac{1}{k}E[X].
$$

Assume $X$ is continuous with PDF $f(\cdot)$. The case for a discrete $X$ is similar.

$$
\begin{aligned}
	E[X] &= \int_0^\infty xf(x)dx \\\\
	&\geq \int_k^\infty xf(x)dx \\\\
	&\geq \int_k^\infty kf(x)dx \quad\cdots\text{ because } x \geq k, \\\\
	&= k \int_k^\infty f(x)dx \\\\
	&= kP(X > k) \\\\
	\frac{1}{k}E[X] &\geq P(X > k)
\end{aligned}
$$

Now we move on to the proof of the theorem. Fix $\epsilon > 0$. Note that

$$
\begin{aligned}
	P(|\hat\theta_n - \theta| > \epsilon) &= P\left( (\hat\theta_n - \theta)^2 > \epsilon^2 \right) \\\\
	&\leq \frac{1}{\epsilon^2} \underbrace{E\left[ (\hat\theta_n - \theta)^2 \right]}\_{MSE \rightarrow 0} \rightarrow 0
\end{aligned}
$$

### Example using the MSE theorem
$Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$. $\hat\theta_n = \max(Y_1, \cdots, Y_n)$. We want to show $\hat\theta_n$ is consistent.

We have the same setup as the [uniform distribution example](#uniform-distribution-example), but this time we want to apply the theorem. Recall that

$$
MSE(\hat\theta_n) = \frac{2\theta^2}{(n+1)(n+2)} \rightarrow 0
$$

as $n \rightarrow \infty$.

### Law of large numbers

The LLN states that if $Y_i$ are i.i.d. with $E[Y_i] = \mu$ and $Var(Y_i) = \sigma^2 < \infty$, then $\bar{Y}_n \xrightarrow{P} \mu$.

**Proof:** $MSE(\bar{Y}_n) = Var(\bar{Y}_n) = \frac{\sigma^2}{n} \rightarrow 0$ as $n \rightarrow \infty$.



{{% alert note %}}

The assumption $Var(Y_i) < \infty$ can be relaxed by assuming only $E[|X_i|] < \infty$.

{{% /alert %}}



## Properties of convergence in probability

Suppose random variables $X_n \xrightarrow{P} x$ and $Y_n \xrightarrow{P} y$.

- $X_n + Y_n \xrightarrow{P} x + y$.
- $X_n Y_n \xrightarrow{P} xy$.
- $X_n / Y_n \xrightarrow{P} x / y$ if $Y_n$ and $y \neq 0$.
- If $g$ is a continuous function, $g(X_n) \xrightarrow{P} g(x)$.

The proofs can be found in more advanced courses. An example of applying them would be $X_n \xrightarrow{P} x$, $Y_n \xrightarrow{P} y$ and $Z_n \xrightarrow{P} z$. We then have

$$
(X_n + Y_n)e^{Z_n} \xrightarrow{P} (x+y)e^z.
$$

Think of this as $X_n' = X_n + Y_n$ and $Y_n' = e^{Z_n}$.

### Example
Suppose $Y_i$ are i.i.d. with $E[Y_i] = \mu$, $Var(Y_i) = \sigma^2 < \infty$ (and $E[Y_i^4] < \infty$). Show that

$$
\hat\sigma_n^2 = \frac{1}{n-1}\sum_{i=1}^n(Y_i - \bar{Y}_n)^2
$$

is unbiased for $\sigma^2$.

The goal is to show $\hat\sigma_n^2$ is consistent for $\sigma^2$. There are several approaches. The straightforward method is to compute $MSE = Var$. This is nasty!

The second (easier) approach is to write

$$
\begin{aligned}
	\hat\sigma_n^2 &= \frac{1}{n-1}\left( \sum_{i=1}^n Y_i^2 - n\bar{Y}_n^2 \right) \\\\
	&= \frac{n}{n-1}\left( \frac{1}{n}\sum_{i=1}^n Y_i^2 - \bar{Y}_n^2 \right) \\
\end{aligned}
$$

As $n \rightarrow \infty$, $\frac{n}{n-1} \rightarrow 1$, $\bar{Y}_n^2 \rightarrow \mu^2$ (LLN), and if we think of $Y_i^2$ as $X_i$,


$$
\frac{1}{n}\sum_{i=1}^n Y_i^2 \xrightarrow{LLN} E[Y_i^2] = \mu^2 + \sigma^2
$$


So we have
$$
\hat\sigma_n^2 \rightarrow 1 \times (\mu^2 + \sigma^2 - \mu^2) = \sigma^2
$$

A note of caution is that the concept of consistency only tells you the convergence eventually. It doesn't tell us how *fast* it's happening. If it's really slow, then it still might not be a good estimator.