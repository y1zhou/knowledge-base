---
title: "Maximum Likelihood Estimator"
slug: "mathematical-statistics-maximum-likelihood-estimator"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Estimation

summary: ""
date: 2020-01-29T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 140

menu:
  maths-stat:
    name: Maximum Likelihood Estimator
    parent: Estimation Under Parametric Models
    weight: 10

---



So far when we talk about estimators, we're using bias, variance etc. to evaluate them. These concepts can be used in general, but if we constrain ourselves within a class of models, we can get much more mathematically useful results.



## Basic concepts

**Def:** Let $\Theta$ be the parameter space (could be high dimensional). $F_\theta$ is a probability  distribution which is uniquely determined by $\theta \in \Theta$. $\\{F_\theta: \theta \in \Theta\\}$ is a `parametric family` of distributions (parameterized by $\theta \in \Theta$).

For example, recall a Poisson distribution with mean $\theta$ is given by PMF

$$
p(x; \theta) = e^{-\theta}\frac{\theta^x}{x!},\quad x = 0, 1, 2, \cdots
$$

If we collect all the PMFs together, $\\{p(x; \theta), \theta \in \Theta\\} = (0, \infty)$ is a parametric family.

Let's consider a negative example. Our model is $Y_i$ are i.i.d. with mean $E[Y_i] = \mu$, which is our parameter.

$$
\\{\text{All distributions with mean }\mu, \mu \in \mathbb{R}\\}
$$

is **not** a parametric family, because given a particular $\mu$, there are multiple distributions having this $\mu$ as a mean.

Another example is $\vec{\theta} = (\mu, \sigma^2)$. Note that we consider $\sigma^2$ as a parameter, not $\sigma$. $\vec{\theta} \in \Theta = (-\infty, \infty) \times (0, \infty)$ [^Cartesian product].

$$
\{N(\mu, \sigma^2): (\mu, \sigma^2) \in \Theta \}
$$

is a parametric family.

**Remark:** it can happen that multiple parameter values correspond to the same distribution (not `identifiable`). For example, $\{N(|\theta|, 1), \theta \in \mathbb{R}\}$. What we cannot have is for one parameter value to map to multiple distributions. We often only consider the identifiable case in this course.

So why is this useful? Given $Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} F_\theta$ where $\theta$ is unknown, we want to find $\theta$ which makes $F_\theta$ fit the data well.

**Def:** statistical analysis based on modeling data with a parametric family is called `parametric statistics`.

## Maximum likelihood estimator

**Def:** Suppose $Y_1, \cdots, Y_n$ are i.i.d. samples that are either

- from a discrete distribution with PMF $p(y; \theta)$, or
- from a continuous distribution with PDF $f(y; \theta)$

where $\theta \in \Theta$ and $\\{p(y, \theta): \theta \in \Theta\\}$ or $\\{f(y, \theta): \theta \in \Theta\\}$ forms a parametric family. The `likelihood function` $L(\theta) = L(\theta; Y_1, \cdots, Y_n)$ is defined as

$$
\begin{gather*}
	\text{discrete}: L(\theta) = p(Y_1; \theta)p(Y_2; \theta)\cdots p(Y_n; \theta) \\\\
	\text{continuous}: L(\theta) = f(Y_1; \theta)f(Y_2; \theta)\cdots f(Y_n; \theta)
\end{gather*}
$$

In other words, a likelihood function is the joint PMF or PDF with the random sample plugged in the variables. Although the likelihood function is essentially a joint PMF/PDF, we are emphasizing on the dependence on the parameter $\theta$.

Given the likelihood function $L(\theta)$, the `maximum likelihood estimator` (MLE) is given by

$$
\begin{equation} \label{eq:MLE}
    \hat\theta = \underset{\theta \in \Theta}{\arg\max}\, L(\theta)
\end{equation}
$$

i.e. find $\theta \in \Theta$ which maximizes $L(\theta)$.

**Remark:** $\hat\theta$ is a statistic that only depends on $Y_1, \cdots, Y_n$. The idea of the MLE is to choose the parameter which makes the sample most probable/likely.

The $\hat\theta$ found in $\eqref{eq:MLE}$ is the same as

$$
\underset{\theta \in \Theta}{\arg\max}\, \ln(L(\theta))
$$

because the $\ln(\cdot)$ function is strictly monotonically increasing. This is called the `log-likelihood function`. This transforms the product in the definition of the likelihood function into a summation, which makes taking the derivative much easier.

$$
\ln(L(\theta)) = \begin{cases}
	\ln(p(Y_1; \theta)) + \cdots + \ln(p(Y_n; \theta)), \text{ or } \\\\
	\ln(f(Y_1; \theta)) + \cdots + \ln(f(Y_n; \theta))
\end{cases}
$$

### Normal distribution example
Suppose $Y_1, \cdots, Y_n \overset{i.i.d.}{\sim}N(\theta, 1)$ and $\theta \in \Theta = \mathbb{R}$. To find the likelihood, we first need the marginal PDF which is given by

$$
f(y; \theta) = \frac{1}{\sqrt{2\pi}} \exp\left\\{ -\frac{(y-\theta)^2}{2} \right\\}
$$

The log-transform of this is

$$
\ln(f(y; \theta)) = -\frac{1}{2}\ln(2\pi) - \frac{(y-\theta)^2}{2}
$$

where the first part is a constant w.r.t. $\theta$. The log-likelihood is

$$
\begin{aligned}
	\ell(\theta) &= \ln(L(\theta)) = \sum_{i=1}^n \ln(f(Y_i; \theta)) \\\\
	&= -\frac{n}{2}\ln(2\pi) - \frac{1}{2}\sum_{i=1}^n (Y_i - \theta)^2
\end{aligned}
$$

Taking the first derivative and setting it to 0, we have
$$
\begin{aligned}
    \frac{d}{d\theta}\ell(\theta) &= 0 - \frac{1}{2}\sum_{i=1}^n 2(Y_i - \theta)(-1) \\\\
	&= \sum_{i=1}^n (Y_i - \theta) = \sum_{i=1}^n Y_i - n\theta = 0
\end{aligned}
$$


which gives $\hat\theta = \bar{Y}_n$. Rigorously speaking, we should then take the second-order derivative to confirm that $\hat\theta$ is the maximizer (get the concave shape). Here this step is omitted.

Note that the maximum likelihood estimator under this model led to the least squares estimator, where we minimized $\sum_{i=1}^n (Y_i - \theta)^2$.

If both parameters in $N(\mu, \sigma^2)$ are unknown, one can take partial derivatives to show that the MLE is given by

$$
(\hat\mu, \hat\sigma^2) = \left( \bar{Y}_n, \frac{1}{n}\sum\_{i=1}^n (Y_i - \bar{Y}_n)^2 \right)
$$

### Bernoulli distribution example
Suppose $Y_i \overset{i.i.d.}{\sim} Bern(\theta)$ and $\theta \in \Theta = [0, 1]$. The marginal PMF is given by

$$
p(y; \theta) = \theta^y (1-\theta)^{1-y}, \quad y = 0, 1
$$

Log-transforming this yields

$$
\ln(p(y; \theta)) = y\ln(\theta) + (1-y)\ln(1-\theta)
$$

The log-likelihood is given by

$$
\begin{aligned}
	\ell(\theta) &= \sum_{i=1}^n \ln(p(Y_i; \theta)) \\\\
	&= \sum_{i=1}^n \bigg(y\ln(\theta) + (1-y)\ln(1-\theta)\bigg) \\\\
	&= \underbrace{X_n\ln(\theta) + (n - X_n)\ln(1-\theta)}_{\text{maximize w.r.t. } \theta} \text{ where } X_n = \sum_{i=1}^n Y_i
\end{aligned}
$$

In the first case where $X_n = 0$, $\ell(\theta) = n\ln(1-\theta)$ is a monotonically decreasing function in $\theta \in [0, 1]$. The maximum is achieved at $\theta = 0$. In this case, the MLE $\hat\theta = 0$.

Similarly in the second case where $X_n = n$, $\ell(\theta) = X_n \ln(\theta)$ and the MLE $\hat\theta = 1$.

When $0 < X_n < n$, we'll have to take the derivative.

$$
\begin{gather*}
	\frac{d}{d\theta} \ell(\theta) = X_n \frac{1}{\theta} - (n - X_n)\frac{1}{1-\theta} = 0 \\\\
	\frac{X_n}{\theta} = \frac{n - X_n}{1 - \theta} \\\\
	X_n - X_n \theta = n\theta - X_n \theta \\\\
	\hat\theta = \frac{X_n}{n} = \bar{Y}_n
\end{gather*}
$$

Going back to the first two cases, we can easily see that $\bar{Y}_n$ also incorporates the two cases, so $\bar{Y}_n$ is the MLE.

### Uniform distribution example
Suppose $Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$ where $\theta \in (0, \infty)$. Find the MLE.

The marginal PDF is

$$
f(y; \theta) = \begin{cases}
	\frac{1}{\theta}, & 0 \leq y \leq \theta \\\\
	0, & \text{otherwise}
\end{cases}
$$

To solve this problem, we first need to introduce the `indicator` $I\{S\}$. It's defined as 

$$
I\{S\} = \begin{cases}
	1 & \text{if } S \text{ holds}, \\\\
	0 & \text{if } S \text{ does not hold}
\end{cases}
$$

where $S$ is a statement. Now we may rewrite the marginal PDF as

$$
f(y; \theta) = \frac{1}{\theta}I\\{y \in [0, \theta]\\}
$$

The likelihood function is given by

$$
\begin{aligned}
	L(\theta) &= f(Y_1; \theta)f(Y_2; \theta) \cdots f(Y_n; \theta) \\\\
	&= \frac{1}{\theta^n}I\\{Y_1 \in [0, \theta]\\}\cdot I\\{Y_2 \in [0, \theta]\\}\cdots I\\{Y_n \in [0, \theta]\\}
\end{aligned}
$$

A useful fact is that

$$
I\\{S_1\\}\cdot I\\{S_2\\} \cdots I\\{S_n\\} = I\\{S_1, \cdots, S_n \text{ all hold}\\}
$$

so the product of the terms reduces to


$$
\begin{aligned}
	&I\\{0 \leq Y_i \leq \theta \text{ for all } i = 1, \cdots n\\} \\\\
	&= I\\{Y_i \leq \theta \text{ for all } i = 1, \cdots n\\} \\\\
	&=  I\\{\max(Y_1, \cdots, Y_n) \leq \theta\\} \\\\
	&=Y_{(n)}
\end{aligned}
$$

Going back to the likelihood function,

$$
\begin{gather*}
	L(\theta) = \frac{1}{\theta^n} I\\{\theta \geq Y_{(n)}\\} \\\\
	\hat\theta = Y_{(n)}
\end{gather*}
$$
because $L(\theta) = 0$ when $\theta < Y_{(n)}$, and is monotonically decreasing after $\theta = Y_{(n)}$.

**Remark:** Often the maximizing problem in the MLE cannot be solved analytically. In this case, we have to go to numerical methods such as gradient descent.

## Reparameterization

Suppose $\\{F_\theta : \theta \in \Theta\\}$ is a parametric family which is identifiable. We have a one-to-one mapping $f: \Theta \rightarrow \Theta'$. If $\hat\theta'$ is the MLE for $\\{F_{f(\theta')}: \theta' \in \Theta'\\}$, then $\hat\theta = f(\hat\theta')$ is the MLE for $\{F_\theta : \theta \in \Theta\}$. MLE is invariant with respect to re-parameterization.

For example, suppose $Y_i \overset{i.i.d.}{\sim} Bern(p)$ where $p \in [0, 1]$. If we're interested in the MLE of $q = 1 - p$, then

$$
MLE(q) = 1 - \hat{p} = \bar{Y}_n
$$

To understand why the MLE is a good estimate, see [Fisher information](https://en.wikipedia.org/wiki/Fisher_information), which can be seen as the curvature of the plot of the log-likelihood.

$$
I(\theta) = -E\left[ \frac{d^2}{d\theta} \ln(p(Y_1; \theta)) \right] \xlongequal{\text{if }\theta \text{ is true}} Var\left( \frac{d}{d\theta} \ln(p(Y_1; \theta)) \right)
$$

Suppose $\theta$ is true. Under some conditions, the MLE $\hat\theta_n$ of $\theta$ satisfies the following theorems:

1. $\hat\theta_n \xrightarrow{P} \theta$ (consistency).
2. $\hat\theta_n$ is approximately distributed as $N \left(\theta, \frac{1}{nI(\theta)} \right)$ when $n$ is large.

The interpretation of the second theorem are

1. $\hat\theta_n$ is asymptotically unbiased.
2. MSE $\approx$ variance $\approx \frac{1}{nI(\theta)}$. View this as two parts: $1/n$ means as the sample size increases, the variance becomes smaller; $1/I(\theta)$ is the curvature, and a large curvature is better.
3. No other consistent estimator can beat this asymptotically. This is called the `Cramér–Rao lower bound `.



[^Cartesian product]: https://en.wikipedia.org/wiki/Cartesian_product