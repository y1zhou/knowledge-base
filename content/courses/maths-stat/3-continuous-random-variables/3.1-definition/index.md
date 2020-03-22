---
title: "Definitions for Continuous Random Variables"
slug: "mathematical-statistics-continuous-rv-definition"
categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Random Variable

summary: "The probability density function, cumulative distribution function, expectation and variance for a continuous random variable."
date: 2019-09-25T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 50

menu:
  maths-stat:
    name: Definitions
    identifier: continuous-definitions
    parent: Continuous Random Variables
    weight: 10
---

In the [previous chapter]({{< ref "/courses/maths-stat/2-discrete-random-variables/2.1-definition/index.md" >}}) we mainly focused on discrete random variables whose set of possible values is either finite or countably infinite. In this chapter, we study random variables whose set of possible values is uncountable. We'll see later on that a lot of the cases we've discussed have analogs in the continuous case.

## Probability density function

The continuous random variable is a random variable with infinite possible outcomes (a subset of the real line). We say $X$ is a `continuous random variable` if there exists a nonnegative function $f$ defined for all $x \in (-\infty, \infty)$, having the property that for any set $B$ of real numbers,
$$
P\\{x \in B\\} = \int_B f(x)dx
$$
Here the function $f$ is called the `probability density function` (PDF). It resembles the probability mass function in the discrete case. The PDF has the following properties:

1. $\int_{-\infty}^\infty f(x)dx = P\\{X \in (-\infty, \infty)\\} = 1$
2. $P\\{a \leq x \leq b\\} = \int_a^b f(x)dx$
3. $P\\{X = a\\} = \int_a^a f(x)dx = 0$

We can also define the `cumulative distribution function` for a continuous random variable:
$$
\begin{aligned}
    F(a) &= P\\{X \leq a\\}, \quad B \in (-\infty, a] \\\\
    &= \int_{-\infty}^a f(x)dx
\end{aligned}
$$

### Properties example

These properties often come in handy when we have unknown quantities in a PDF. Suppose $X$ is a continuous random variable with probability density function
$$
f(x) = \begin{cases}
    C\left( 4x - 2x^2 \right), & 0 < x < 2 \\\\
    0, & \text{otherwise}
\end{cases}
$$
and we'd like to find $C$ as well as the probability $P\\{X > 1\\}$.
$$
\begin{aligned}
    1 &= \int_{-\infty}^\infty f(x)dx \\\\
    &= \int_{-\infty}^0 f(x)dx + \int_0^2 f(x)dx + \int_2^\infty f(x)dx \\\\
    &= C\int_0^2 (4x - 2x^2)dx \\\\
    &= C\left(\int_0^2 4xdx - \int_0^2 2x^2dx \right) \\\\
    &= C \left( 2x^2\bigg|_0^2 - \frac{2}{3}x^3 \bigg|_0^2 \right) \\\\
    &= C\left( 8 - \frac{16}{3} \right)
\end{aligned}
$$
Now we can easily find $C = \frac{3}{8}$. With the PDF given, it's trivial to find the CDF:


$$
\begin{aligned}
    P\{X > 1\} &= \int_1^\infty f(x)dx \\\\
    &= \int_1^2 \frac{3}{8}\left(4x - 2x^2\right)dx \\\\
    &= \frac{3}{8}\left( 2x^2\bigg|_1^2 - \frac{2}{3}x^3 \bigg|_1^2 \right) \\\\
    &= \frac{3}{8}\left( 8-2 - \left(\frac{16}{3} - \frac{2}{3}\right) \right) = \frac{1}{2}
\end{aligned}
$$


### Lifetime example

Suppose $X$, the lifetime of an item, is a continuous random variable  with a density function
$$
f(x) = \begin{cases}
    \lambda e^{-x/100}, & x \geq 0 \\\\
    0, & x < 0
\end{cases}
$$
What is the probability that the item functions between 50 and 150 days?

In mathematical terms, we want to calculate $P\\{50 \leq X \leq 150\\}$. We first need to find the value of $\lambda$.


$$
\begin{aligned}
    1 &= \int_{\infty}^\infty f(x)dx \\\\
    &= \int_0^\infty \lambda e^{-\frac{x}{100}}dx \\\\
    &= \lambda \int_0^\infty e^{-\frac{x}{100}}dx \\\\
    &= -100\lambda \int_0^\infty e^{-\frac{x}{100}} d\frac{-x}{100} \\\\
    &= -100\lambda e^{-\frac{x}{100}} \bigg|_0^\infty \\\\
    &= -100\lambda(0 - 1) = 100\lambda
\end{aligned}
$$


> Recall that $\int e^xdx = e^x$, and $d(ax) = a \cdot dx$ because the derivative is a linear function.

With $\lambda = \frac{1}{100}$, we can calculate


$$
\begin{aligned}
    P\{50 \leq X \leq 150\} &= \int_{50}^{150} \frac{1}{100}e^{-\frac{x}{100}}dx \\\\
    &= \frac{-100}{100}\int_{50}^{150}e^{-\frac{x}{100}} d\frac{-x}{100} \\\\
    &= -\left( e^{-\frac{x}{100}}\bigg|_{50}^{150} \right) \\\\
    &= e^{-\frac{1}{2}} - e^{-\frac{3}{2}} \approx 0.383
\end{aligned}
$$


## Expectation and variance

Earlier we've defined the expectation for discrete random variables. If $X$ is a continuous random variable with probability density function $f(x)$, we have
$$
f(x)dx \approx P\\{x \leq X \leq x + dx\\}
$$


so it's easy to find the analog for the `expectation` of $X$ to be


$$
E[X] = \int_{-\infty}^\infty xf(x)dx
$$


Similarly, the expected value of a real-valued function of $X$ is


$$
E[g(x)] = \int_{-\infty}^\infty g(x)f(x)dx
$$


which can be used to derive the `variance` of $X$
$$
Var(X) = E[X^2] - E[X]^2
$$


Suppose $X$ is a continuous random variable with density function
$$
f(x) = \begin{cases}
    2x, & 0 \leq x \leq 1 \\\\
    0, & \text{otherwise}
\end{cases}
$$

$$
\begin{aligned}
    E[X] &= \int_{-\infty}^\infty xf(x)dx \\\\
    &= \int_0^1 x \cdot 2x dx \\\\
    &= \frac{2}{3}x^3 \bigg|_0^1 = \frac{2}{3} \\\\
    E[X^2] &= \int_{-\infty}^\infty x^2 f(x)dx \\\\
    &= \int_0^1 2x^3 dx \\\\
    &= \frac{1}{2}x^4 \bigg|_0^1 = \frac{1}{2} \\\\
    Var(X) &= \frac{1}{2} - \left(\frac{2}{3}\right)^2 = \frac{1}{18}
\end{aligned}
$$


Next, we introduce some commonly seen continuous probability distributions.