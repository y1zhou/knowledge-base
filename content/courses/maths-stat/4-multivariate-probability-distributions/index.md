---
title: "Multivariate Probability Distributions"
slug: "mathematical-statistics-multivariate-probability-distributions"
categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Random Variable

summary: "Joint probability distributions of two or more random variables defined on the same sample space. Also covers independence, conditional expectation and total expectation."
date: 2019-11-06T09:57:16-04:00
toc: true
type: docs  # Do not modify.
weight: 70

menu:
  maths-stat:
    name: Multivariate Probability Distributions
    # identifier: YourPageID
    # parent: YourParentID
    weight: 40
---

In the previous two chapters, we focused on studying a single random variable. However, we can define more than one random variable on the same sample space.

Let's consider the following motivating example. Suppose we roll two six-sided fair dice. The sample space contains 36 possible sample points. Define


$$
\begin{aligned}
    &X_1 \text{: the number of die 1} \\\\
    &X_2 \text{: the number of die 2} \\\\
    &X_3 \text{: the sum of the two dice} \\\\
    &X_4 \text{: the product of the two dice}
\end{aligned}
$$
Sometimes, we may want to assign a probability to an event that involves two or more random variables. For example,


$$
\begin{aligned}
    P\\{X_1 = 1, X_2 = 2\\} &= \frac{|E|}{|S|} = \frac{|E_1 \cap E_2|}{|S|} \\\\
    &= \frac{|(1, 2)|}{|S|} = \frac{1}{36}
\end{aligned}
$$
This is called the `joint probability` of $X_1$ and $X_2$.



## Two discrete random variables

**Def:** Let $X_1$ and $X_2$ be two discrete random variables defined on the same sample space. The `joint probability mass function` for $X_1$ and $X_2$ is given by

$$
p(x_1, x_2) = P\\{X_1 = x_1, X_2 = x_2\\}
$$
The joint probability mass function has the following properties:

1. $p(x_1, x_2) \geq 0$ for all $x_1$ and $x_2$.
2. $\sum\limits_{X_2} \sum\limits_{X_1} p(x_1, x_2) = 1$.

Given the joint PMF, we can also define the `joint cumulative distribution function` (also called the joint probability distribution function) as
$$
\begin{aligned}
    F(a, b) &= P\\{X_1 \leq a, X_2 \leq b\\} \\\\
    &= \sum_{X_1 \leq a, X_2 \leq b} p(x_1, x_2) \\\\
    &= \sum_{X_1 \leq a} \sum_{X_2 \leq b} p(x_1, x_2)
\end{aligned}
$$


To distinguish with the joint probability mass function and the joint probability distribution function, we call the mass and distribution functions of *a single random variable as* `marginal` mass and distribution functions.

Suppose we have two discrete random variables $X$ and $Y$,
$$
\begin{aligned}
    p_X(x) &= P\\{X = x\\} \\\\
    &= P\left\\{ (X = x) \cap S \right\\} \\\\
    &= P\left\\{ \\{X = x\\} \cap \\{Y \leq \infty\\} \right\\} \\\\
    &= \sum_{y: p(x, y) > 0} p(x, y)
\end{aligned}
$$
The marginal probability mass function of $X$ can be obtained by summing the joint probability mass function over possible values of $Y$. The marginal distribution function can be derived from the joint distribution function:
$$
F_X(a) = P\\{X \leq a\\} = P\\{X \leq a, Y < \infty\\} = F_{X, Y}(a, \infty)
$$
Suppose we have a box that contains 3 red balls and 4 blue balls. Let
$$
\begin{aligned}
    &X \text{: number of red balls,} \\\\
    &Y \text{: number of blue balls}
\end{aligned}
$$
If we randomly draw 3 balls out of the 7 balls, find the joint probability mass function $p(x, y)$.

We can list out all possible configurations of this random experiment:



| $_{Y = j} \backslash^{X = i}$ | 0    | 1    | 2    | 3    |
| ----------------------------- | ---- | ---- | ---- | ---- |
| 0                             | 0    | 0    | 0    |      |
| 1                             | 0    | 0    |      | 0    |
| 2                             | 0    |      | 0    | 0    |
| 3                             |      | 0    | 0    | 0    |



The blanks in the table are the possible events. We have $X + Y = 3$, so the possible outcomes are
$$
\begin{gather*}
    \\{ X = 3, Y = 0 \\} \\\\
    \\{ X = 2, Y = 1 \\} \\\\
    \\{ X = 1, Y = 2 \\} \\\\
    \\{ X = 0, Y = 3 \\}
\end{gather*}
$$
We can find their probabilities


$$
\begin{aligned}
    P\\{ X = 3, Y = 0 \\} &= p(3, 0) = \frac{|\\{X = 3, Y = 0\\}|}{|S|} \\\\
    &= \frac{\binom{3}{3}\binom{4}{0}}{\binom{7}{3}} = \frac{1}{35} \\\\
    P\\{ X = 2, Y = 1 \\} &= p(2, 1) = \frac{\binom{3}{2} \binom{4}{1}}{\binom{7}{3}} = \frac{12}{35}
\end{aligned}
$$
Similarly we can find $p(1, 2) = \frac{18}{35}$ and $p(0, 3) = \frac{4}{35}$. We can also find the marginal probability of $X$ at, say $X = 2$:


$$
\begin{aligned}
    p_X(2) &= P\\{X = 2\\} \\\\
    &= P\\{X = 2, Y \leq \infty\\} \\\\
    &= \sum_{y: p(x, y) > 0} P\\{X = 2, Y = y\\} \\\\
    &= p(2, 1) = \frac{12}{35}
\end{aligned}
$$


## Two continuous random variables

Let $X$ and $Y$ be two continuous random variables. We say $X$ and $Y$ are `jointly continuous` if there exists a function $f(x, y) \geq 0$, such that


$$
P\\{(x, y) \in B\\} = \iint\limits_{(x, y) \in B} f(x, y)dxdy
$$


The function $f(x, y)$ is called the `joint probability density function` of $X$ and $Y$. We can decompose the event
$$
\\{(X, Y) \in B\\} \text{ as } \\{X \in B_X, Y \in B_Y\\}
$$
So we can also rewrite the joint density function into


$$
P\\{X \in B_X, Y \in B_Y\\} = \int\limits_{Y \in B_Y} \int\limits_{X \in B_X} f(x, y)dxdy
$$
The `joint cumulative distribution function` can be defined using integrals of the joint density functions:


$$
\begin{aligned}
    F(a, b) &= P\\{X \leq a, Y \leq b\\} \\\\
    &= \int\limits_{-\infty}^b \int\limits_{-\infty}^a f(x, y)dxdy \\\\
    &= \int\limits_{-\infty}^a \int\limits_{-\infty}^b f(x, y)dxdy
\end{aligned}
$$
Also, we can find the relationship between the joint and `marginal density functions`:


$$
\begin{aligned}
    f_X(x) = \int\limits_{-\infty}^\infty f(x, y)dy \\\\
    f_Y(y) = \int\limits_{-\infty}^\infty f(x, y)dx
\end{aligned}
$$
Now, suppose the joint density function of continuous random variables $X$ and $Y$ is given by


$$
f(x, y) = \begin{cases}
    2e^{-x}e^{-2y}, & x > 0, y > 0 \\\\
    0, & \text{otherwise}
\end{cases}
$$
and we want to compute

1. $P\\{X > 1, Y < 1\\}$,
2. $P\\{X < Y\\}$, and
3. $P\\{X < a\\}$.

For case 1,
$$
\begin{aligned}
    P\\{X > 1, Y < 1\\} &= \int\limits_{-\infty}^1 \int\limits_1^\infty f(x, y)dxdy \\\\
    &= \int\limits_0^1 \int\limits_1^\infty 2e^{-x}e^{-2y}dxdy \\\\
    &= \int\limits_0^1 2e^{-2y}dy \int\limits_1^\infty e^{-x}dx \\\\
    &= \int\limits_0^1 2e^{-2y}dy \left( -e^{-x}\bigg|_1^\infty \right) \\\\
    &= e^{-1} \int\limits_0^1 2e^{-2y}dy \\\\
    &= e^{-1} \left( -e^{-2y}\bigg|_0^1 \right) \\\\
    &= e^{-1}\left( 1 - e^{-2} \right)
\end{aligned}
$$
For case 2,
$$
\begin{aligned}
    P\\{X < Y\\} &= P\\{ -\infty < y < \infty, x < y \\} \\\\
    &= \int\limits_{-\infty}^\infty \int\limits_{x < y} f(x, y)dxdy \\\\
    &= \int\limits_0^\infty \int\limits_0^y 2e^{-x}e^{-2y}dxdy \\\\
    &= \int\limits_0^\infty 2e^{-2y}dy \int\limits_0^y e^{-x}dx \\\\
    &= \int\limits_0^\infty 2e^{-2y}dy \left( 1 - e^{-y} \right) \\\\
    &= \int\limits_0^\infty 2e^{-2y}dy - 2\int\limits_0^\infty e^{-3y}dy \\\\
    &= e^{-2y} \bigg|_\infty^0 + \frac{2}{3} e^{-3y} \bigg|_0^\infty \\\\
    &= 1 - \frac{2}{3} = \frac{1}{3}
\end{aligned}
$$
And for case 3,
$$
\begin{aligned}
    P\\{X < a\\} &= \int\limits_{-\infty}^a f_X(x)dx \\\\
    &= \int\limits_{-\infty}^a \int\limits_{-\infty}^\infty f(x, y)dydx \\\\
    &= \int\limits_0^a \int\limits_0^\infty 2e^{-x}e^{-2y}dydx \\\\
    &= \int\limits_0^a e^{-x}dx \underbrace{\int\limits_0^\infty 2e^{-2y}dy}_{Exp(2)} \\\\
    &= \int\limits_0^a e^{-x}dx = F_X(a) \\\\
    &= e^{-x}\bigg|_a^0 \\\\
    &= 1 - e^{-a}
\end{aligned}
$$


## Conditional probability distribution

If we have events $E$ and $F$ on a sample space $S$, the [conditional probability]({{< ref "/courses/maths-stat/1-probability/1.2-conditional-probability/index.md" >}}) of $E$ given $F$ is defined as
$$
P(E \mid F) = \frac{P(EF)}{P(F)}
$$


We can easily extend this idea to two random variables. Suppose we now have random variables $X$ and $Y$, and we define events $E = \\{X = x\\}$ and $F = \\{Y = y\\}$, then


$$
P\\{E \mid F\\} = P\\{ X = x \mid Y = y \\} = \frac{P\\{X = x, Y = y\\}}{P\\{Y = y\\}}
$$


### The discrete case

**Def:** If $X$ and $Y$ are two discrete random variables, we can define the `conditional probability mass function` of $X$ given $Y$ as

$$
P_{X \mid Y}(x \mid y) = P\\{ X = x \mid Y = y \\} = \frac{p(x, y)}{P_Y(y)}
$$


Similarly, the `conditional distribution function` of $X$ given $Y = y$ can be defined as


$$
F_{X \mid Y}(a \mid y) = P\\{X \leq a \mid Y = y\\} = \sum\limits_{X \leq a} P_{X \mid Y}(x \mid y)
$$


We can also check if the conditional mass function satisfies the properties of a probability mass function:


$$
\begin{aligned}
    P_{X \mid Y}(x \mid y) &= g(x), \quad \text{treating }y \text{ as constant} \\\\
    &= \frac{p(x, y)}{P_Y(y)} \\\\
    \sum\limits_X \frac{p(x, y)}{P_Y(y)} &= \frac{1}{P_Y(y)} \sum\limits_X p(x, y) \\\\
    &= \frac{1}{P_Y(y)} P_Y(y) = 1
\end{aligned}
$$


As an example, suppose that $p(x, y)$ is given by
$$
\begin{gather*}
    p(0, 0) = 0.4, && p(0, 1) = 0.2, \\\\
    p(1, 0) = 0.1, && p(1, 1) = 0.3
\end{gather*}
$$
We want to calculate the conditional probability mass function of $X$ given that $Y = 1$. Using the notations above,


$$
P_{X \mid Y}(X \mid 1) = \frac{p(x, 1)}{P_Y(1)}
$$
We have
$$
p(x, 1) = \begin{cases}
    p(0, 1) = 0.2, & X = 0 \\\\
    p(1, 1) = 0.3, & X = 1
\end{cases}
$$
and the marginal probability function of $Y$
$$
\begin{aligned}    P_Y(1) &= \sum\limits_X p(x, 1) \\    &= p(0, 1) + p(1, 1) \\    &= 0.2 + 0.3 = 0.5\end{aligned}
$$
The conditional PMF is given by


$$
P_{X \mid Y}(X \mid 1) = \begin{cases}
    \frac{0.2}{0.5} = \frac{2}{5}, & X = 0 \\\\
    \frac{0.3}{0.5} = \frac{3}{5}, & X = 1
\end{cases}
$$


### The continuous case

**Def:** If $X$ and $Y$ are two continuous random variables, we can define the `conditional joint density function` of $X$ given $Y = y$, for all values of $y$ such that $f_Y(y) > 0$, by

$$
f_{X \mid Y}(x \mid y) = \frac{f(x, y)}{f_Y(y)}
$$


The conditional distribution function of $X$ given $Y = y$ can be calculated as


$$
\begin{aligned}
    F_{X \mid Y}(a \mid y) &= P\{X \leq a \mid Y = y\} \\\\
    &= \int\limits_{-\infty}^a f_{X \mid Y}(x \mid y)dx
\end{aligned}
$$
Again we'll explain through an example. Suppose the joint density function of continuous random variables $X$ and $Y$ is given by
$$
f(x, y) = \begin{cases}
    \frac{12}{5}x(2 - x - y), & 0< x < 1, 0 < y < 1 \\\\
    0, & \text{otherwise}
\end{cases}
$$
We want to compute the conditional density of $X$ given $Y = y$ when $0 <y < 1$.


$$
\begin{aligned}
    f_{X \mid Y}(x \mid y) &= \frac{f(x, y)}{f_Y(y)}, \quad 0<y<1 \\\\
    f_Y(y) &= \int\limits_{-\infty}^\infty f(x, y)dx \\\\
    &= \int\limits_0^1 \frac{12}{5}x(2-x-y)dx \\\\
    &= \left( \frac{6}{5}x^2(2-y) - \frac{4}{5}x^3 \right) \bigg|_0^1 \\\\
    &= -\frac{4}{5} + \frac{6}{5}(2-y) \\\\
    &= \frac{12}{5}\left( \frac{2-y}{2} - \frac{1}{3} \right) \\\\
    &= \frac{12}{5}\left( \frac{2}{3} - \frac{y}{2} \right)
\end{aligned}
$$

$$
f_{X \mid Y}(x \mid y) = \begin{cases}
    \frac{6x(2-x-y)}{4 - 3y}, & \text{when } 0<x<1 \\\\
    0, & \text{otherwise}
\end{cases} \quad 0 < y < 1
$$


## Independent random variables

The concept of [independence between events]({{< ref "/courses/maths-stat/1-probability/1.2-conditional-probability/index.md#independence-of-events" >}}) can also be extended to that between two random variables. If $X$ and $Y$ are two random variables, we say $X$ and $Y$ are `independent` with each other, i.e. $X \perp Y$, if for all $A$ and $B$, the events $E_A = \\{X \in A\\}$ and $E_B = \\{Y \in B\\}$ are independent with each other: 


$$
P\\{X \in A, Y \in B\\} = P\\{X \in A\\}P\\{Y \in B\\} = P\\{E_A\\}P\\{E_B\\}
$$
If we take $E_A = \\{X \leq a\\}$ and $E_B = \\{Y \leq b\\}$ for a pair of real numbers $a$ and $b$,
$$
P\\{X \leq a, Y \leq b\\} = P\\{X \leq a\\} P\\{Y \leq b\\}
$$
Hence, in terms of the joint distribution function of $X$ and $Y$,


$$
F_{X, Y}(a, b) = F_X(a)F_Y(b) \quad \forall a, b \in \mathbb{R}
$$


If $X$ and $Y$ are **discrete random variables** and $X \perp Y$, the condition of independence is equivalent to


$$
\begin{gather*}
    p(x, y) = P_X(x)P_y(y) \\\\
    P_{X \mid Y}(x \mid y) = P_X(x) \\\\
    P_{Y \mid X}(y \mid x) = P_Y(y)
\end{gather*}
$$


In the **continuous case**,


$$
X \perp Y \Leftrightarrow \begin{cases}
    f(x, y) = f_X(x)f_Y(y) \\\\
    f_{X \mid Y}(x \mid y) = f_X(x) \\\\
    f_{Y \mid X}(y \mid x) = f_Y(y)
\end{cases}
$$


So loosely speaking, $X$ and $Y$ are independent if knowing the value of one doesn't change the distribution of another.

As an example, suppose a couple decided to date at a restaurant for dinner. Suppose each of them independently arrives at a time uniformly distributed between 6 p.m. and 7 p.m. Find the probability that the first to arrive has to wait longer than 10 minutes.

Let $X$ and $Y$ denote the number of minutes past 6 that the man and woman arrives, respectively. We have


$$
X \sim Unif(0, 60), \quad Y \sim Unif(0, 60)
$$
and $X \perp Y$. Define events


$$
E_1 = \\{Y > X + 10\\} \text{ and } E_2 = \\{X > Y + 10\\}
$$
What we want to find is $P\{E_1 \cup E_2\}$, that is the man ($E_1$) or the woman ($E_2$) waits for more than 10 minutes. Since $E_1$ and $E_2$ are disjoint, $P\{E_1 \cup E_2\} = P(E_1) + P(E_2)$. By symmetry,


$$
P(E_1) = P(E_2) \Rightarrow P\{E_1 \cup E_2\} = 2P(E_1)
$$
We may express $E_1$ as


$$
E_1 = \\{Y > X + 10\\} = \\{ 10<Y \leq 60, 0 \leq X < Y - 10 \\}
$$
To find its probability,


$$
\begin{aligned}
    P(E_1) &= \int\limits_{10}^{60} \int\limits_0^{y-10} f(x, y)dxdy \\\\
    &= \int\limits_{10}^{60} \int\limits_0^{y-10}f_X(x) f_Y(y)dxdy
\end{aligned}
$$
where the equality holds because of independence.


$$
\begin{aligned}
    P(E_1) &= \int\limits_{10}^{60} \int\limits_0^{y-10} \frac{1}{60 - 0} \cdot \frac{1}{60}dxdy \\\\
    &= \frac{1}{3600}\int\limits_{10}^{60} dy x \bigg|_0^{y-10} \\\\
    &= \frac{1}{3600} \int\limits_{10}^{60} (y-10)dy \\\\
    &= \frac{1}{3600} \left( \frac{y^2}{2} - 10y \right) \bigg|_{10}^{60} \\\\
    &= \frac{25}{72}
\end{aligned}
$$
So the final answer is $\frac{25}{36}$.



## Expectation

Recall that for a random variable $X$, its expectation is given by


$$
\begin{aligned}
    E[X] &= \sum\limits_X xP_X(x) \quad \text{for a discrete r.v.} \\\\
    &= \int\limits_{-\infty}^\infty xf_X(x)dx \quad \text{for a continuous r.v.}
\end{aligned}
$$


and both can be considered as special cases of $E[g(X)]$:


$$
\begin{aligned}
    E[g(X)] &= \sum\limits_X g(x)P_X(x) \quad \text{for a discrete r.v.} \\\\
    &= \int\limits_{-\infty}^\infty g(x)f_X(x)dx \quad \text{for a continuous r.v.}
\end{aligned}
$$

**Def:** If $X$ and $Y$ are two **discrete random variables**, we can extend the definition using the joint probability mass function:

$$
E[g(X, Y)] = \sum\limits_X \sum\limits_Y g(x, y) p(x, y)
$$


and similarly using the joint probability density function for **the continuous case**:


$$
E[g(X, Y)] = \int\limits_{-\infty}^\infty \int\limits_{-\infty}^\infty g(x, y)f(x, y)dxdy
$$


Suppose $X$ and $Y$ have the joint density function


$$
f(x, y) = \begin{cases}
    4xy, & 0 \leq x \leq 1, 0 \leq y \leq 1 \\\\
    0, & \text{otherwise}
\end{cases}
$$


The expectation is given by


$$
\begin{aligned}
    E[XY] &= \int\limits_{-\infty}^\infty \int\limits_{-\infty}^\infty \underbrace{xy}_{g(X, Y)} f(x, y) dxdy \\\\
    &= \int\limits_0^1 \int\limits_0^1 xy \cdot 4xy dxdy \\\\
    &= \int\limits_0^1 4y^2dy \int\limits_0^1 x^2dx \\\\
    &= \int\limits_0^1 4y^2dy \left(\frac{x^3}{3}\right) \bigg|_0^1 \\\\
    &= \frac{4}{3} \int\limits_0^1 y^2 dy = \frac{4}{9}
\end{aligned}
$$


### Properties

1. Similar to the property of the expectation of a single random variable,

$$
E[ag(x, y) + bh(x, y)]  = aE[g(x, y)] + bE[h(x, y)]
$$



2. If $g(x, y) > h(x, y)$, then $E[g(x, y)] > E[h(x, y)]$.

3. If $X \perp y$ and $g(X, Y) = g_1(X)g_2(Y)$,

$$
E[g(x, y)] = E[g_1(X)g_2(Y)] = E[g_1(X)] E[g_2(Y)]
$$

If we can decompose the joint density function $f(x, y)$, we can show that it's the production of the two marginal probability functions. A special case of (3) is $E[XY] = E[X] E[Y]$.

These properties could make a lot calculations much easier. Suppose $X$ and $Y$ are two standard normal random variables and are independent with each other and we want to calculate $E\left[(2X + Y)^2\right]$.

With the properties shown above, we no longer need to find $f(x, y) = f_X(x) f_Y(y)$.


$$
\begin{aligned}
    E\left[(2X + Y)^2\right] &= E\left[ 4X^2 + 4XY + Y^2 \right] \\\\
    &= 4E\left[X^2\right] + 4E[XY] + E\left[ Y^2 \right] \\\\
    &= 4E\left[X^2\right] + 4E[X]E[Y] + E\left[ Y^2 \right] \\\\
    E\left[X^2\right] &= Var(X) + E[X]^2 = 1 \\\\
    E\left[(2X + Y)^2\right] &= 4 + 0 + 1 = 5
\end{aligned}
$$


## Covariance

As shown in the example above, we have defined the expected value and variance for a single random variable. For two or more random variables, we introduce a new quantity to give us the information about the relationship between the random variables.

**Def:** The `covariance` of $X$ and $Y$, denoted as $Cov(X, Y)$, is defined as
$$
Cov(X, Y) = E\bigg[ (X - E[X])(Y - E[Y]) \bigg]
$$
Upon expanding the right hand side of the definition above,
$$
\begin{aligned}
    Cov(X, Y) &= E[XY - XE[Y] - E[X]Y + E[X]E[Y]] \\\\
    &= E[XY] - E[Y]E[X] - E[X]E[Y] + E[X]E[Y] \\\\
    &= E[XY] - E[X]E[Y]
\end{aligned}
$$


### Properties

1. The order doesn't matter:

$$
Cov(X, Y) = Cov(Y, X)
$$



2. The relationship between the covariance and the variance of a single random variable:

$$
Cov(X, X) = Var(X)
$$



3. Similar to the property of the variance,

$$
Cov(aX + b, cY + d) = ac \cdot Cov(X, Y)
$$



4. For $X_1, \cdots, X_n$ and $Y_1, \cdots, Y_n$,

$$
Cov\left( \sum_{i=1}^n X_i, \sum_{j=1}^n Y_j \right) = \sum_{j=1}^n \sum_{i=1}^n Cov(X_i, Y_j)
$$



5. For independent random variables,

$$
X \perp Y \Rightarrow Cov(X, Y) = 0
$$

Note that the arrow doesn't go both ways! Having zero covariance (`uncorrelated`) doesn't imply independence. An example is $X \sim Unif(-1, 1)$ and $Y = X^2$.

As an example, suppose $X$ and $Y$ are two continuous random variables with joint density function


$$
f(x, y) = \begin{cases}
    3x, & 0 \leq y \leq x \leq 1 \\\\
    0, & \text{otherwise}
\end{cases}
$$


and we want to find $Cov(X, Y)$.

We can view the joint density function as $g_1(X) = 3x$, $y \leq x \leq 1$ and $g_2(Y) = 1$, $0 \leq y \leq x$. To find $Cov(X, Y)$, we need to find $E[XY]$, $E[X]$ and $E[Y]$.


$$
\begin{aligned}
    E[XY] &= \int\limits_{-\infty}^\infty \int\limits_{-\infty}^\infty xy f(x, y)dxdy \\\\
    &= \int\limits_0^1 \int\limits_y^1 xy 3x dxdy \\\\
    &= \int\limits_0^1 3ydy \int\limits_y^1 x^2 dx \\\\
    &= \int\limits_0^1 3ydy \left(\frac{x^3}{3}\right) \bigg|_y^1 \\\\
    &= \int\limits_0^1 3ydy \left( \frac{1}{3} - \frac{y^3}{3} \right) \\\\
    &= \int\limits_0^1 ydy - \int\limits_0^1 y^4dy \\\\
    &= \left( \frac{y^2}{2} - \frac{y^5}{5} \right) \bigg|_0^1 \\\\
    &= \frac{3}{10}
\end{aligned}
$$
To calculate $E[X]$, we first need to find the marginal probability function of $X$.


$$
\begin{aligned}
    f_X(x) &= \int\limits_{-\infty}^\infty f(x, y)dy \\\\
    &= \int\limits_0^x 3xdy \\\\
    &= 3x^2, \quad 0 \leq x \leq 1 \\\\
    E[X] &= \int\limits_0^1 x \cdot 3x^2 dx \\\\
    &= \frac{3}{4}x^4 \bigg|_0^1 = \frac{3}{4}
\end{aligned}
$$


And similarly for $E[Y]$,


$$
\begin{aligned}
    f_Y(y) &= \int\limits_{-\infty}^\infty f(x, y) dx \\\\
    &= \int\limits_y^1 3xdx = \frac{3}{2}x^2 \bigg|_y^1 \\\\
    &= \frac{3}{2}(1 - y^2), \quad 0 \leq y \leq 1 \\\\
    E[Y] &= \int\limits_0^1 y \cdot \frac{3}{2}(1 - y^2)dy \\\\
    &= \frac{3}{2} \left( \frac{y^2}{2} - \frac{y^4}{4} \right) \bigg|_0^1 \\\\
    &= \frac{3}{4} - \frac{3}{8} = \frac{3}{8}
\end{aligned}
$$
Now we can show the covariance between $X$ and $Y$ to be
$$
Cov(X, Y) = E[XY] - E[X]E[Y] = \frac{3}{10} - \frac{3}{4} \cdot \frac{3}{8} = \frac{3}{160}
$$


## Conditional expectation

Recall that we have the conditional probability function defined as


$$
P\left( \underbrace{X \in B_X}_E \mid \underbrace{Y \in B_Y}_F \right) = \frac{P(X \in B_X, Y \in B_Y)}{P(Y \in B_Y)}
$$


[In the discrete case](#the-discrete-case), we can define the `conditional expectation` of $X$ given $Y = y$ for all values of $y$ such that $p_Y(y) > 0$ by


$$
E[X \mid Y = y] = \sum_X xP_{X \mid Y}(x \mid y)
$$


In the [continuous case](#the-continuous-case),


$$
E[X \mid Y = y] = \int\limits_{-\infty}^\infty x f_{X \mid Y}(x \mid y) dx = \int\limits_{-\infty}^\infty x \frac{f(x, y)}{f_Y(y)}dx
$$


### Discrete example

Suppose $X$ and $Y$ are independent binomial random variables with parameters $n$ and $p$. Find the conditional expected value of $X$ given that $X + Y = m$.


$$
\begin{aligned}
    E[X \mid X + Y = m] &= \sum_X x P_{X \mid (X + Y)}(x \mid m) \\\\
    P_{X \mid (X + Y)}(x \mid m) &= P \\{ X = x \mid (X + Y) = m \\} \\\\
    &= \frac{P \\{ X = x, X + Y = m \\}}{P \\{ X + Y = m \\}}
\end{aligned}
$$
Since $X = 0, 1, 2, \cdots, n$ and $X + Y = m$, if $m < n$, the values $X$ can take would be limited to $X = 0, 1, 2, \cdots, m$. For $k \leq \min(n, m)$,


$$
\begin{aligned}
    P_{X \mid X + Y}(x \mid m) &= \frac{P \\{ X = k, X + Y = m \\}}{P \\{ X + Y = m \\}}  \\\\
    &= \frac{P\\{ X = k, Y = m-k \\}}{P\\{X + Y = m\\}} \\\\
    &= \frac{P\\{ X = k\\} P\\{ Y = m-k \\}}{P\\{X + Y = m\\}} \\\\
    P\\{X = k\\} &= p_X(k) = \binom{n}{k}p^k(1-p)^{n-k} \\\\
    P\\{Y = m-k\\} &= p_Y(m-k) = \binom{n}{m-k}p^{m-k}(1-p)^{n-m+k} \\\\
    P\\{X + Y = m\\} &= P\\{Z = m\\}, \quad Z \sim Bin(2n, p) \\\\
    &= \binom{2n}{m}p^m(1-p)^{2n-m} \\\\
    P_{X \mid X + Y}(x \mid m) &= \frac{\binom{n}{k}p^k(1-p)^{n-k}\binom{n}{m-k}p^{m-k}(1-p)^{n-m+k}}{\binom{2n}{m}p^m(1-p)^{2n-m}} \\\\
    &= \frac{\binom{n}{k}\binom{n}{m-k}}{\binom{2n}{m}}
\end{aligned}
$$


Here $Z \sim Bin(2n, p)$ because $X$ and $Y$ are independent, and we can consider $X + Y$ as $2n$ independent trails, each with a success probability $p$.

Now we can plug the $P_{X \mid X + Y}(x \mid m)$ term back into the conditional expectation equation:


$$
\begin{equation} \label{eq:conditional-expectation}
    E[X \mid X + Y = m] = \sum_{k=0}^m k \frac{\binom{n}{k}\binom{n}{m-k}}{\binom{2n}{m}}
\end{equation}
$$
Note that another way of writing out the probability mass function of $P\\{X + Y = m\\}$ is
$$
\begin{aligned}
    P\\{X + Y = m\\} &= \sum_{k=0}^m P\\{X = k\\}P\\{Y = m-k\\} \\\\
    &= \sum_{k=0}^m \binom{n}{k}p^k(1-p)^{n-k} \binom{n}{m-k}p^{m-k}(1-p)^{n-m+k} \\\\
    &= \sum_{k=0}^m \binom{n}{k} \binom{n}{m-k} p^m(1-p)^{2n-m} \\\\
    &= \binom{2n}{m}p^m(1-p)^{2n-m} \\\\
    \Rightarrow \binom{2n}{m} &= \sum_{k=0}^m \binom{n}{k} \binom{n}{m-k}
\end{aligned}
$$




Now we need to rewrite the binomial terms in Equation $\eqref{eq:conditional-expectation}$ to remove $k$:


$$
\begin{aligned}
    k \binom{n}{k} &= k \frac{n!}{k!(n-k)!} = \frac{n!}{(k-1)!(n-k)!} \\\\
    &= \frac{n(n-1)!}{(k-1)![(n-1) - (k-1)]!} = n \binom{n-1}{k-1} \\\\
    \binom{n}{m-k} &= \binom{n}{(m-1)-(k-1)} \\\\
    \binom{2n}{m} &= \frac{(2n)!}{m!(2n-m)!} \\\\
    &= \frac{2n(2n-1)!}{m(m-1)![(2n-1)-(m-1)]!} \\\\
    &= \frac{2n}{m}\binom{2n-1}{m-1}
\end{aligned}
$$


Going back to Equation $\eqref{eq:conditional-expectation}$,


$$
\begin{aligned}
    E[X \mid X + Y = m] &= \sum_{k=0}^m \frac{n\binom{n-1}{k-1}\binom{n}{(m-1)-(k-1)}}{\frac{2n}{m}\binom{2n-1}{m-1}} \\\\
    &= \frac{m}{2} \frac{1}{\binom{2n-1}{m-1}}\sum_k \binom{n-1}{k-1}\binom{n}{(m-1)-(k-1)} \\\\
    &= \frac{m}{2}
\end{aligned}
$$


### Continuous example

$X$ and $Y$ are continuous random variables with joint density function


$$
f(x, y) = \begin{cases}
    \frac{e^{-\frac{x}{y}}e^{-y}}{y}, & 0 < x < \infty, 0 < y < \infty \\\\
    0, & \text{otherwise}
\end{cases}
$$
Find $E[X \mid Y = y]$ where $0 < y < \infty$.

Note that in the equation of the conditional density


$$
f_{X \mid Y}(x \mid y) = \frac{f(x, y)}{f_Y(y)}
$$
the joint density function is given, so we only need to figure out the marginal density function of $Y$.


$$
\begin{aligned}
    f_Y(y) &= \int\limits_{-\infty}^\infty f(x, y)dx \\\\
    &= \int\limits_0^\infty \frac{e^{-\frac{x}{y}}e^{-y}}{y}dx \\\\
    &= \frac{e^{-y}}{y} \int\limits_0^\infty e^{-\frac{x}{y}}dx \\\\
    &= -e^{-y} \int\limits_0^\infty e^{-\frac{x}{y}}d\left(-\frac{x}{y}\right) \\\\
    &= -e^{-y} \left(e^{-\frac{x}{y}}\right)\bigg|_0^\infty \\\\
    &= e^{-y}, \quad 0 < y < \infty
\end{aligned}
$$


So the conditional density function is
$$
\begin{aligned}
    f_{X \mid Y}(x \mid y) &= \frac{e^{-\frac{x}{y}}e^{-y}y^{-1}}{e^{-y}} \\\\
    &= \frac{1}{y}e^{-\frac{x}{y}}, \quad 0 < x < \infty, 0 < y < \infty
\end{aligned}
$$
which follows an exponential density function with parameter $y$. Finally we can find the conditional expectation


$$
E[X \mid Y = y] = \int\limits_0^\infty x \cdot \frac{1}{y}e^{-\frac{x}{y}}dx = y
$$


because for $Z \sim Exp(y)$, $f(Z) = \frac{1}{y} e^{-\frac{Z}{y}}$ and $E[Z] = y$.



### Properties

1. Summation comes outside:

$$
E\left[\sum_{i=1}^n X_i \mid Y = y\right] = \sum_{i=1}^n E[X_i \mid Y = y]
$$



2. 

$$
E[aX \mid Y = y] = aE[X \mid Y = y]
$$



3. Combining 1 and 2,

$$
E \left[ \sum_{i=1}^n a_iX_i \mid Y = y \right] = \sum_{i=1}^n a_i E\left[ X_i \mid Y = y \right]
$$

4. Generalizing 2,

$$
E[g(X) \mid Y = y] = \begin{cases}
    \sum_X g(X)P_{X \mid Y}(x \mid y), & \text{discrete r.v.s} \\\\
    \int_{-\infty}^\infty  g(X) f_{X \mid Y}(x \mid y), & \text{continuous r.v.s}
\end{cases}
$$

### Making predictions

A very common problem in statistics is to make predictions. Suppose the value of a random variable $X$ is observed, and an attempt is made to predict the value of another random variable $Y$ based on this observed value of $X$.

Let $g(X)$ denote the predictor. Clearly, we would want a function $g(\cdot)$ such that $g(X)$ tends to be close to $Y$. A popular criterion for closeness is the error sum of squares. The model is given by


$$
\arg \min\limits_{g(\cdot)} E\left[ (Y - g(X))^2 \right]
$$
We can think of this as $Y = g(X) + \epsilon$ where $\epsilon$ is the estimation error. The equation above is minimizing the variance of $\epsilon$. We can show that under this criterion, the best possible predictor of $Y$ is $g(X) = E[Y \mid X]$.

**Theorem:** $E\left[ (Y - g(X))^2 \right] \geq E\left[ (Y - E[Y \mid X = x])^2 \right]$. The proof is given below.


$$
\begin{aligned}
    E\left[ (Y - g(X))^2 \mid X = x \right] &= E\left[ \bigg(\underbrace{Y - E[Y \mid X]}_{a} + \underbrace{E[Y \mid X] - g(X)}_b\bigg)^2 \mid X \right] \\\\
    &= E\left[ \left(a^2 + b^2 + 2ab\right) \mid X \right] \\\\
    &= E[a^2 \mid X] + E[b^2 \mid X] + 2E[ab \mid X]
\end{aligned}
$$

$$
\begin{aligned}
    E[ab \mid X] &= E\bigg[ (Y - E[Y \mid X])(\underbrace{E[Y \mid X]}_{h(X)} - g(X)) \mid X \bigg] \\\\
    &= \Big( E[Y \mid X] - g(X) \Big) E\Big[Y - E[Y \mid X] \mid X \Big] \\\\
    &= \Big( E[Y \mid X] - g(X) \Big) \Big( E[Y \mid X]- E[Y \mid X] \Big) = 0
\end{aligned}
$$

$$
\begin{gather*}
    E\left[ (Y - g(X))^2 \mid X = x \right] = E[a^2 \mid X] + E[b^2 \mid X] \\\\
    E\left[ (Y - g(X))^2 \right] = E\left[ (Y - E[Y \mid X])^2\right] + \underbrace{E[b^2]}_{b \geq 0} \\\\
    E\left[ (Y - g(X))^2 \right] \geq E\left[ (Y - E[Y \mid X])^2\right]
\end{gather*}
$$


Now suppose the height of a man is $X$ inches and the height of his son is $Y$ inches Assume $Y \sim N(X + 1, 4)$. If we know a father is 6 feet tall, what is the expected value of the height of his son?

The model can be written as


$$
Y = X + 1 + \epsilon
$$


where $\epsilon \sim N(0, 4)$ and is independent of $X$. The best prediction of $Y$ given $X = x$ is $E[Y \mid X = x]$.


$$
\begin{aligned}
    E[Y \mid X = 72] &= E[X + 1 + \epsilon \mid X = 72] \\\\
    &= 73 + E[\epsilon \mid X = 72] \\\\
    &= 73 + E[\epsilon] = 73
\end{aligned}
$$


## Total expectation

For any value of $y$, we can think of the conditional expectation as a real-valued function of $Y$,


$$
E[X \mid Y = y] = h(Y)
$$
and $h(Y) = H$ can be considered as another random variable.


$$
\begin{equation} \label{eq:total-expectation}
    E[H] = \underbrace{E\bigg[ \underbrace{E[X \mid Y = y]}_{\text{about } Y} \bigg]}_{\text{about } X} = E[X]
\end{equation}
$$


The equation $E[X] = E[E[X \mid Y = y]]$ is called `total expectation`. If $X$ and $Y$ are two **discrete random variables**,


$$
E[X] = \sum_y \underbrace{E[X \mid Y = y]}_{h(y)} p_Y(y)
$$


In the **continuous case**,


$$
E[X] = \int\limits_{-\infty}^\infty \underbrace{E[X \mid Y = y]}_{h(y)} f_Y(y) dy
$$


The proof for the discrete case of Equation $\eqref{eq:total-expectation}$ is given below. The proof for the continuous case is similar.


$$
\begin{aligned}
    \sum_y E[X \mid Y = y]p_Y(y) &= \sum_y \left( \sum_X xP_{X \mid Y}(x \mid y) \right) p_Y(y) \\\\
    &= \sum_y \left( \sum_X x \underbrace{P_{X \mid Y}(x \mid y)p_Y(y)}_{P_{X, Y}(x, y)} \right) \\\\
    &= \sum_y \sum_X x P_{X, Y}(x, y) \\\\
    &= \sum_X x \sum_y P_{X, Y}(x, y) \\\\
    &= \sum_X xp_X(x) = E[X]
\end{aligned}
$$


**Example:** Suppose Joe is in a room with doors A, B and C. Door A leads outside and it takes 3 hours to go through the tunnel behind it. Doors B and C both lead back to the room, and the tunnels take 5 and 7 hours to finish, respectively. The three doors are shuffled each time one door is chosen, so a random guess is needed each time. What is the expected value for the number of hours to leave the room?

Let $X$ be the number of hours to escape. The possible values of $X$ goes from 3 to infinity (if Joe's really unlucky). Directly calculating $E[X]= \sum_X xp_X(x)$ is therefore not straightforward, so we define


$$
\begin{aligned}
    Y = \\{\text{the door initially chosen}\\} = \quad 
    &1 && 2 && 3 \\\\
    &A && B && C
\end{aligned}
$$
and calculate $E[X]$ using total expectation.


$$
\begin{aligned}
    E[X] &= E\bigg[ E[X \mid Y = y] \bigg] \\\\
    &= \sum_y E[X \mid Y = y]p_Y(y) \\\\
    &= E[X \mid Y = 1]p_Y(1) + E[X \mid Y = 2]p_Y(2) + E[X \mid Y = 3]p_Y(3) \\\\
    &= \frac{1}{3} \bigg( E[X \mid Y = 1] + E[X \mid Y = 2] + E[X \mid Y = 3] \bigg) \\\\
    &= \frac{1}{3} \bigg(3 + (E[X] + 5) + (E[X] + 7) \bigg)
\end{aligned}
$$


which gives us $E[X] = 15$. To understand the conditions, we explain $E[X \mid Y = 2]$ as an example. If Joe chose the second door, he would spend 5 hours in the tunnel and then return to the room. Since the doors are shuffled, the problem is the same as before. Thus, his expected additional time until leaving the room is just $E[X]$, and $E[X \mid Y = 2] = 5 + E[X]$.