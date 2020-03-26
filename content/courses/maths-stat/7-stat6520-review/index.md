---
title: "Brief Review Before STAT 6520"
slug: "mathematical-statistics-brief-review-before-6520"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics

summary: "A brief review of probability theory "
date: 2020-01-08T09:57:16-04:00
toc: true
type: docs  # Do not modify.
weight: 100

menu:
  maths-stat:
    name: Brief Review
    # identifier: YourPageID
    # parent: YourParentID
    weight: 70
---



{{% alert note %}}

This chapter gives a brief review of STAT 6510 before we move on to STAT 6520. It's poorly written now and needs major rework.

{{% /alert %}}

Why do we learn probability theory first? In practice, it's usually impossible or undesirable to establish a deterministic relation in observed data. There is randomness involved in the data. To deal with this randomness, we have to relax the deterministic relation by some **unexplainable** randomness. In statistical analysis, we assume that the `data generating mechanism` (model) involves randomness.



For the first part of this course, we'll focus on a very simple type of model: i.i.d. sample $Y_1, \cdots, Y_n$ from the same (population) distribution, i.e. $Y_i \overset{i.i.d.}{\sim} F$.



>  Notation: $X \sim F$ means random variable $X$ follows distribution $F$, e.g. $X \sim N(0, 1)$.



Suppose we toss the same coin a number of times and record the outcomes. Let $Y_i$ denote the outcome of the $i$-th coin toss, and
$$
Y_i = \begin{cases}
    0, & \text{Tail} \\\\
    1, & \text{Head}
\end{cases}
$$


**Model:** $Y_i \overset{i.i.d.}{\sim} Bern(p)$ where $0 < p < 1$. The probability mass function of each $Y_i$ is
$$
f(y) = P(Y_i = y) = \begin{cases}
    1-p, & y = 0 \\\\
    p, & y = 1
\end{cases}
$$




Now we can ask more questions such as

- what are $E[Y_i]$ and $Var(Y_i)$? Since $Y_i$ is a discrete random variable, we have


$$
\begin{gather*}
    E[Y_i] = \sum_y yf(y) =0(1-p) + 1 \cdot p = p \\\\
    Var(Y_i) = E\left[ (Y_i - E[Y_i])^2 \right] = E[Y_i^2] - E[Y_i]^2 \\\\
    E[Y_i^2] = \sum_y y^2f(y) = 1^2 \times p = p \\\\
    Var(Y_i) = p - p^2 = p(1-p)
\end{gather*}
$$


Taking a step further, define $X_n = Y_1 + Y_2 + \cdots + Y_n$, which is the total number of heads after $n$ tosses.

- What is the distribution of $X_n$?

  $X_n \sim Bin(n, p)$. The probability mass function of $X_n$ is


$$
g(x) = P(X_n = x) = \binom{n}{x}p^x(1-p)^{n-x}, \quad x = 0, 1, 2, \cdots, n
$$


Its expectation is


$$
\begin{aligned}
    E[X_n] &= E[Y_1 + \cdots + Y_n] \\\\
    &= E[Y_1] + \cdots + E[Y_n] \\\\
    &= np
\end{aligned}
$$




And its variance is


$$
\begin{aligned}
    Var(X_n) &= Var(Y_1 + \cdots + Y_n) \\\\
    &= Var(Y_1) + \cdots + Var(Y_n) \quad\cdots\text{independence} \\\\
    &= np(1-p)
\end{aligned}
$$


- Now, what is the behavior of $\frac{X_n}{n}$, the frequency of heads, as $n \rightarrow \infty$?

  Intuitively, the answer is $p$ (by the law of large numbers). 

- What continuous distribution does $X_n$ look like when $n$ is large?

   The distribution of $X_n$ will become very close to a normal distribution with parameters $np$ and $np(1-p)$. This is the central limit theorem.



All questions above are treated in probability theory, where we study the properties of the (random) data generated from the given model. What is the goal of `statistics`? We still have a model $Y_i \overset{i.i.d.}{\sim} Bern(p)$, but we assume that $p$ is **unknown** and we want to make inference about $p$. 



## Law of large numbers and the central limit theorem

The [law of large numbers](https://en.wikipedia.org/wiki/Law_of_large_numbers) (LLN) states that if we have i.i.d. random variables $Y_i$, as the sample size $n \rightarrow \infty$, the sample mean approaches the population mean:

$$
\frac{1}{n}\sum_{i=1}^n Y_i \rightarrow E[Y_i]
$$

The [central limit theorem](https://en.wikipedia.org/wiki/Central_limit_theorem) (CLT) states that if we sum the i.i.d. random variables $Y_i$ to get $X_n = Y_1 + \cdots + Y_n$, the distribution of $X_n$ is approximately $N(\mu, \sigma^2)$ where $\mu = E[X_n] = nE[Y_i]$ and $\sigma^2 = Var(X_n) = nVar(Y_i)$.

In the example of i.i.d. Bernoulli random variables above, by the `LLN` as $n$ increases, the distribution of $X_n$ shifts its center to the right; by the CLT the distribution of $\frac{X_n}{n}$ concentrates to the center $p$ as $n$ increases (higher peak).



## Probability theory and statistics

In probability theory, we're given a model and we generate random data/sample. We then study the properties of the data. In statistics, we're given a model with an **unspecified** element $p$. After the random data/sample is generated, we loop back to the model and try to make inference to the unknown parameter $p$. The two fundamental concepts in statistical inference are `estimation` and `hypothesis testing`.

1. Provide an estimation of $p$, e.g. $\hat{p} = \frac{X_n}{n}$. How good is this estimation?

2. Is the hypothesis $p = 0.5$ supported by the data? How do we make the decision?



{{% alert warning %}}

All the statistical theories we will develop in this course are based on certain model assumptions. However, any model assumption is an approximation of the reality. As George Box stated, all models are wrong but some are useful.

{{% /alert %}}