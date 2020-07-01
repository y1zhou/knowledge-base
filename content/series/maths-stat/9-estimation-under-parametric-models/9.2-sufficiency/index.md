---
title: "Sufficiency"
slug: "mathematical-statistics-sufficiency"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Estimation

summary: "Introducing sufficient statistics for the inference of parameters. The factorization theorem comes in handy!"
date: 2020-01-30T10:46:18-04:00
toc: true
weight: 92
---



{{% alert warning %}}

This is one of the hardest topics in this course.

{{% /alert %}}



We have a single sample $Y \sim N(0, \sigma^2)$ where $\sigma^2$ is unknown. Suppose the sign of $Y$ is lost, so only $|Y|$ is retained. Does this affect the guess of $\sigma^2$?

Not really, because our distribution is centered at $0$, so the sign of the observations doesn't really affect how far it's away from the mean.

## Definition

Suppose $Y_i \overset{i.i.d.}{\sim} Bern(p)$, $p \in [0, 1]$. If only $T = Y_1 + \cdots + Y_n$ is known, does it affect the estimation of $p$?

The answer is again no, because the estimator we use is

$$
\bar{Y}_n = \frac{1}{n}\sum\_{i=1}^n Y_i = \frac{1}{n}T
$$

When we say $T$ is known, we're actually referring to the concept of `conditioning` in probability theory. Consider the conditional distribution

$$
\begin{gather*}
	P(Y_1 = y_1, Y_2 = y_2, \cdots, Y_n = y_n \mid T = t) = \frac{P(A \cap B)}{P(B)} \\\\
	= \frac{P(Y_1 = y_1, \cdots, Y_n = y_n, T = t)}{P(T = t)}
\end{gather*}
$$

where each $y_i = 0$ or $1$. If we know the value of $t$, we actually know how many $1$'s there are in $y_1, \cdots, y_n$.

$$
\begin{aligned}
	&\cdots = \frac{P(Y_1 = y_1, \cdots, Y_n = y_n) I\left\\{ \sum_{i=1}^n y_i = t \right\\}}{P(T = t)} \\\\
	&= \frac{P(Y_1 = y_1)^nI\left\\{ \sum_{i=1}^n y_i = t \right\\}}{P(T = t)} \\\\
	&= \frac{p^{\sum_{i=1}^n Y_i}(1-p)^{n-\sum_{i=1}^n Y_i}I\left\\{ \sum_{i=1}^n y_i = t \right\\}}{P(T = t)} \\\\
	&= \frac{p^{\sum_{i=1}^n Y_i}(1-p)^{n-\sum_{i=1}^n Y_i}I\left\\{ \sum_{i=1}^n y_i = t \right\\}}{\binom{n}{t}p^t(1-p)^{n-t}} \\\\
	&= \begin{cases}
		\frac{p^t(1-p)^{n-t}}{\binom{n}{t}p^t(1-p)^{n-t}} = \frac{1}{\binom{n}{t}}, & \sum_{i=1}^n y_i = t, \\\\
		0, & \text{otherwise}
	\end{cases} \\\\
	&= \frac{1}{\binom{n}{t}} I\left\\{\sum_{i=1}^n y_i = t\right\\}
\end{aligned}
$$

**Remark:** If $T$ is given, then the joint distribution of $(Y_1, \cdots, Y_n)$ will no longer depend on $p$.



**Def:** Suppose $Y_1, \cdots, Y_n$ is a sample from a population distribution with unknown parameter $\theta$. A statistic $T$ is said to be `sufficient` for $\theta$ if the conditional distribution of $(Y_1, \cdots, Y_n)$ given $T$ does **not** depend on $\theta$. The interpretation is that knowing $T$ is sufficient for inferring $\theta$. After conditioning on $T$, the sample no longer reflect information about $\theta$.

Suppose statistic $T$ is sufficient for $\theta$. Let $f$ be an one-to-one known function, then $f(T)$ is still sufficient for $\theta$. For example, in the example above we know that $T = \sum_{i=1}^n Y_i$ is sufficient for $p$. With $f(t) = t/n$, we can show that $\bar{Y}_n$ is also sufficient for $p$.

## Factorization theorem

The statistic $T = f(Y_1, \cdots, Y_n)$ is sufficient for $\theta$ if and only if

$$
L(\theta; Y_1, \cdots, Y_n) = g(T, \theta) \times h(Y_1, \cdots, Y_n)
$$

where $g(T, \theta)$ involves both $T$ and $\theta$, whereas $h(\cdot)$ has no $\theta$ involved.

**Remark:** this factorization is not unique! Often $h$ can be chosen as a constant $1$.

Take the $Y_i \overset{i.i.d.}{\sim} Bern(p)$ example again. We have $T = \sum_{i=1}^n Y_i$.

$$
\begin{aligned}
	L(p; Y_1, \cdots, Y_n) &= p^{\sum_{i=1}^n Y_i}(1-p)^{n-\sum_{i=1}^n Y_i} \\\\
	&= \underbrace{p^T(1-p)^{n-T}}_{g} \\\\
	&= g \times 1
\end{aligned}
$$

So by the factorization theorem, $T$ is sufficient for $p$.

### Normal distribution example

$Y \sim N(0, \sigma^2)$ and our sample size $n = 1$. We want to find if $|Y|$ is sufficient.

$$
\begin{aligned}
	L(\sigma^2; Y) &= \frac{1}{\sigma\sqrt{2\pi}}\exp\left( -\frac{Y^2}{2\sigma^2} \right) \\\\
	&= \frac{1}{\sigma\sqrt{2\pi}}\exp\left( -\frac{|Y|^2}{2\sigma^2} \right) \\\\
	&= g(|Y|, \sigma^2) \times 1
\end{aligned}
$$

We can also check the sufficiency by directly using the definition. To find the distribution of $Y$ conditioning on $|Y|$, see that $(Y \mid |Y|)$ can only take two values.

$$
P(Y = y \mid |Y| = t) = \begin{cases}
	\frac{1}{2}, & y = t \\\\
	\frac{1}{2}, & y = -t
\end{cases}
$$

There's no $\sigma^2$ in this distribution.

Going back to the factorization theorem, when $Y_i$ is discrete, let $p(t)$ denote the marginal PMF of $T$. Assuming the factorization holds,

$$
\begin{aligned}
	p(y_1, \cdots, y_n \mid t) &= P(Y_1 = y_1, \cdots, Y_n = y_n \mid T = t) \\\\
	&= \frac{P(Y_1 = y_1, \cdots, Y_n = y_n, T = t)}{p(t)} \\\\
	&= \frac{g(t, \theta) h(y_1, \cdots, y_n)I\{ f(y_1, \cdots, y_n) = t \}}{p(t)} \\\\
	&= g(t, \theta) \cdot h^*(y_1, \cdots, y_n; t)
\end{aligned}
$$


$$
\begin{gather*}
  \sum_{y_1, \cdots, y_n} p(y_1, \cdots, y_n \mid t) = 1 \Leftrightarrow \sum_{y_1, \cdots, y_n} g(t, \theta) h^*(y_1, \cdots, y_n; t) = 1 \\\\
	g(t, \theta)\sum_{y_1, \cdots, y_n} h^*(y_1, \cdots, y_n; t) = 1 \\\\
	g(t, \theta) = \frac{1}{\sum_{y_1, \cdots, y_n} h^*(y_1, \cdots, y_n; t)}
\end{gather*}
$$


The $g(t, \theta)$ function does not depend on $\theta$, therefore $P(Y_1 = y_1, \cdots, Y_n = y_n \mid T = t)$ does not depend on $\theta$.

### Chi square distribution example

Suppose $Y_i$ are i.i.d. samples from PDF

$$
f(y; \nu) = \frac{1}{2^{\frac{\nu}{2}}\Gamma(\frac{\nu}{2})}y^{\frac{\nu}{2} - 1}e^{-\frac{y}{2}} I\\{y > 0\\}
$$

which is the $\chi^2$ distribution with $\nu$ degrees of freedom, and $\nu = 1, 2, 3, \cdots$ is unknown. The likelihood function is

$$
L(\nu; Y_1, \cdots, Y_n) = \underbrace{\frac{1}{2^{\frac{n\nu}{2}}\Gamma(\frac{\nu}{2})^n}(Y_1\cdots Y_n)^{\frac{\nu}{2} - 1}}\_{g} \underbrace{e^{-\frac{Y_1 + \cdots + Y_n}{2}}}\_{h}
$$

Thus $T = Y_1 Y_2 \cdots Y_n$ is the sufficient statistic for $\nu$.

### Example with two parameters

{{% alert note %}}
In the factorization, both $T$ and $\theta$ can be vectors, which means the theorem can be applied when there are multiple parameters.
{{% /alert %}}

$Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} Unif(a, b)$ where $a < b$. In this case $\vec\theta = (a, b)$. Find a sufficient statistic for $\theta$.

The likelihood function is given by

$$
\begin{aligned}
	L(a, b; Y_1, \cdots, Y_n) &= \left(\frac{1}{b-a}\right)^n I\left\\{ a \leq \min(Y_1, \cdots, Y_n), b \geq \max(Y_1, \cdots, Y_n) \right\\} \\\\
	&= \left(\frac{1}{b-a}\right)^n I\left\\{ a \leq Y_{(1)}, b \geq Y_{(n)} \right\\}
\end{aligned}
$$

See [example in previous page]({{< ref "/series/maths-stat/9-estimation-under-parametric-models/9.1-maximum-likelihood-estimator/index.md#uniform-distribution-example" >}}) for details. The entire RHS would be the $g(\cdot)$ part as they all contain $a$ and $b$. So $T = (Y_{(1)}, Y_{(n)})$ is sufficient for $(a, b)$.

## Proposition

The MLE is always a function of (or based on) a sufficient statistic, because the MLE automatically explores full information in the sample.

By the factorization theorem, if $T$ is sufficient, the likelihood function can be written as

$$
L(\theta) = g(T, \theta)h(Y_1, \cdots, Y_n)
$$

In the MLE, maximizing $L(\theta)$ is the same as

$$
\arg\max_{\theta \in \Theta} g(T, \theta)
$$

which depends on $T$ only.
