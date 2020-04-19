---
title: "Optimal Tests"
slug: "mathematical-statistics-optimal-tests"
categories:
    - Mathematical Statistics
tags:
    - Mathematical Statistics
    - Statistics

summary: ""
date: 2020-04-14T18:26:34-04:00

toc: true
type: docs # Do not modify.
weight: 210

menu:
    maths-stat:
        name: Optimal Tests
        parent: Hypothesis Testing
        weight: 40
---

So far we've learnt various aspects about hypothesis tests. We've seen some examples and discussed ways of constructing tests from CIs and concepts like p-values. Now we move on to the theoretical aspect.

Is there a "best" test among all the possible tests that we can have? To understand the notion of optimality, recall that we have constrained the type I error to be $\leq$ the level $\alpha$. Our goal was to minimize the type II error, or equivalently to maximize the power. The optimality of a test is just some constrained optimization.

## Motivating example

Let's consider the following example. Suppose we have a single sample $Y \sim Unif(0, \theta)$. The hypotheses are

$$
H_0: \theta = 1 \quad vs. \quad \theta > 1
$$

Compare the rejection rules

1. $Y > 1 - \alpha$, and
2. $Y < \alpha$.

Intuitively, the first one is a better choice since it's trying to reflect the hypotheses and reject $H_0$ when the sample is large. However, both rejection rules have level $\alpha$ since

$$
\begin{gather*}
  P_1(Y > 1 - \alpha) = 1 - (1 - \alpha) = \alpha \\\\
  P_1(Y < \alpha) = \alpha
\end{gather*}
$$

using the CDF of the uniform distribution. If we compare the power, since $\frac{Y}{\theta} \sim Unif(0, 1)$, for $\theta > 1$,

$$
\begin{gather*}
    pw_1(\theta) = P_\theta(Y > 1 - \alpha) = P_\theta \left(\frac{Y}{\theta} > \frac{1-\alpha}{\theta} \right) = 1 - \frac{1 - \alpha}{\theta} \\\\
    pw_2(\theta) = P_\theta(Y < \alpha) = P_\theta \left(\frac{Y}{\theta} < \frac{\alpha}{\theta} \right) = \frac{\alpha}{\theta}
\end{gather*}
$$

So for all $\theta > 1$,

$$
pw_1(\theta) - pw_2(\theta) = 1 - \frac{1 - \alpha}{\theta} - \frac{\alpha}{\theta} = 1 - \frac{1}{\theta} > 0
$$

Hence the first rejection rule is preferred. If we think about the problem in general, the challenges are

1. The power $pw(\theta)$ is a function of $\theta \in \Theta_\alpha$. There are different $\theta$ values involved.
2. Optimization is with respect to all possible decision rules.

## Neyman-Pearson lemma

If we eliminate the first challenge above and assume that both hypotheses are simple, then the following theorem can be applied to this simplest scenario.

**Assumptions:** We have hypotheses $H_0$: $\theta = \theta_0$ vs. $H_a: \theta = \theta_a$. The significance level is $\alpha$, and $L(\theta)$ is the likelihood function where $\theta \in \\{\theta_0, \theta_a\\}$.

**Conclusion:** The test (decision rule) which maximizes $pw(\theta_a)$ is given by

$$
T = \frac{L(\theta_0)}{L(\theta_a)}, \quad RR = \\{t: t < k_\alpha\\}
$$

where $k_\alpha$ is chosen to make[^neyman-pearson-level] $P_{\theta_0} (T < k_\alpha) = \alpha$. In other words, we reject $H_0$ when the likelihood at $\theta_0$ is relatively small compared with the likelihood at $\theta_a$.

### Proof

Here we only consider the continuous case. Suppose sample $Y = \\{Y_1, \cdots, Y_n\\} \in \mathbb{R}^n$. Any rejection rule can be formulated using a region $D \subset \mathbb{R}^n$ where we reject $H_0$ if $\\{Y_1, \cdots, Y_n\\} \in D$. For a general test statistic $T = t(Y_1, \cdots, Y_n)$, we have

$$
D = \\{(y_1, \cdots, y_n): t(y_1, \cdots, y_n) \in RR\\}
$$

Suppose $D_1, D_2 \subset \mathbb{R}^n$ and

-   Test 1 is $D_1$, the test given in the theorem.
-   Test 2 is $D_2$, an arbitrary test with level $\leq \alpha$.

Our goal is to show powers $pw_1(\theta_a) \geq pw_2(\theta_a)$. Let

$$
\begin{gathered}
    S_1 = \\{y \in \mathbb{R}^n: 1_{D_1}(y) > 1_{D_2}(y)\\} \\\\
    S_1 = \\{y \in \mathbb{R}^n: 1_{D_1}(y) < 1_{D_2}(y)\\} \\\\
\end{gathered}
$$

Namely, in $S_1$ test 1 rejects $H_0$ but test 2 doesn't; in $S_2$ test 1 doesn't reject $H_0$ but test 2 rejects $H_0$.

Let $L_0 = L(\theta_0)$ and $L_a = L(\theta_a)$. Note that

$$
\begin{aligned}
    &\quad \int_{\mathbb{R}^n} \left[ 1_{D_1}(y) - 1_{D_2}(y) \right] \left[ k_\alpha L_a(y) - L_0(y) \right] dy \\\\
    &= \int_{S_1} (1_{D_1} - 1_{D_2})(k_\alpha L_a - L_0)dy + \int_{S_2} (1_{D_1} - 1_{D_2})(k_\alpha L_a - L_0)dy
\end{aligned}
$$

By the rule of test 1, in $S_1$ we have $k_\alpha L_a > L_0$, which means the first term above is $\geq 0$. Similarly in $S_2$, we have $k_\alpha L_a \leq L_0$, so the second term is also $\geq 0$. By rearranging terms, we have

$$
\begin{equation} \label{eq:neyman-pearson}
    k_\alpha \int_{\mathbb{R}_n} (1_{D_1} - 1_{D_2})L_\alpha dy \geq \int_{\mathbb{R}_n} (1_{D_1} - 1_{D_2})L_0 dy
\end{equation}
$$

Note that

$$
\int_{\mathbb{R}_n} 1_{D_i} L_0(y) dy = \int_{D_i}L_0(y)dy = P_{\theta_0}(Y \in D_i)
$$

is the type I error for test $i$ where $i = 1, 2$. A similar relation holds if $L_0$ is replaced by $L_a$ which leads to powers. So Equation $\eqref{eq:neyman-pearson}$ translates to

$$
k_\alpha \left(pw_1(\theta_a) - pw_2(\theta_a) \right) \geq \alpha - \alpha = 0
$$

### Example

Consider $Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} N(\theta, 1)$ and we want to test

$$
H_0: \theta = \theta_0 \quad vs. \quad H_a: \theta = \theta_a, \quad \theta_0 < \theta_a
$$

The likelihood is given by

$$
L(\theta) = \frac{1}{(2\pi)^\frac{n}{2}} \exp \left \\{ -\frac{\sum_{i=1}^n (Y_i - \theta)^2}{2} \right\\}
$$

Therefore,

$$
\begin{aligned}
    T &= \frac{L(\theta_0)}{L(\theta_a)} = \exp \left\\{ -\frac{\sum_{i=1}^n (Y_i - \theta_0)^2 - \sum_{i=1}^n (Y_i - \theta_a)^2}{2} \right\\} \\\\
    &= \exp \left\\{ (\theta_0 - \theta_a)\sum_{i=1}^n Y_i - \frac{n}{2}\left(\theta_0^2 - \theta_a^2 \right) \right\\}
\end{aligned}
$$

and we reject $H_0$ if $T < k_\alpha$ for some suitable $k_\alpha$. We may observe that {{<hl>}} $T$ is a decreasing function of $\bar{Y}\_n$ since we assumed $\theta_0 < \theta_a$. Hence, we reject $H_0$ if $\bar{Y}\_n > c_\alpha$ for some suitable $c_\alpha$. {{</hl>}}

To determine $c_\alpha$, set

$$
\begin{aligned}
    P_{\theta_0}(\bar{Y}_n > c_\alpha) &= P_{\theta_0} \left( \frac{1}{\sqrt{n}}\sum_{i=1}^n (Y_i - \theta_0) > \sqrt{n}(c_\alpha - \theta_0) \right) \\\\
&= 1 - \Phi\left( \sqrt{n}(c_\alpha - \theta_0) \right) = \alpha
\end{aligned}
$$

which leads to

$$
\sqrt{n}(c_\alpha - \theta_0) = z_\alpha \Rightarrow c_\alpha = \frac{z_\alpha}{\sqrt{n}} + \theta_0
$$

since $\bar{Y}_n$ follows $N\left(\theta, \frac{1}{n} \right)$.

The rejection rule can be equivalently formulated using the log-likelihood, where we reject $H_0$ if

$$
\ln L(\theta_0) - \ln L(\theta_a) < u_\alpha = \ln k_\alpha
$$

which often simplifies calculations.

## Composite alternative

The theorem tells us about the most powerful test in the situation of simple $H_0$ vs. simple alternative $H_a$. In general, we will deal with composite $H_0$ and $H_a$. How do we formulate optimality in this more general setting?

**Def:** Suppose we are given a sample and want to test

$$
H_0: \theta \in \Theta_0 \quad vs. \quad H_a: \theta \in \Theta_a
$$

We say a test is the `uniformly most powerful` (UMP) one if for any $\theta \in \Theta_a$, the power $pw(\theta)$ is the largest among all possible choices of decision rules[^hidden-constraint].

**Corollary:** In the case of a simple $H_0: \theta = \theta_0$ and a composite $H_a: \theta \in \Theta_a$, if the decision rule specified in the Neyman-Pearson lemma

$$
\text{reject } H_0 \text{ when } \frac{L(\theta_0)}{L(\theta)} < k_\alpha, \quad \theta \in \Theta_a
$$

**does not** depend on $\theta$, then it yields a UMP test.

To understand this, consider [the example above](#example) where we replace the alternative with $H_a: \theta > \theta_0$, and still use the likelihood ratio statistic. The derivation shows that no matter which $\theta_a > \theta_0$ we focus on, we get the same rejection rule:

$$
\bar{Y}\_n > c_\alpha = \frac{z_\alpha}{\sqrt{n}} + \theta_0
$$

which doesn't have $\theta_a$ involved, hence the test is UMP. This is actually how the Neyman-Pearson lemma typically gets applied. Unfortunately, the applicability of this corollary is very limited. Whenever the alternative is two-sided, often the UMP does not exist.

If the alternative is $H_a: \theta \neq \theta_0$. We know the most powerful (MP) test for any $\theta > \theta_0$ is to reject $H_0$ when

$$
\bar{Y}\_n > \frac{z_\alpha}{\sqrt{n}} + \theta_0
$$

Similarly we can derive the MP test for any $\theta < \theta_0$ is to reject $H_0$ when

$$
\bar{Y}\_n < -\frac{z_\alpha}{\sqrt{n}} + \theta_0
$$

The two rejection rules are mutually exclusive, so there's no UMP test for both sides.

To have a well-defined optimality problem for two-sided tests, one often introduces an additional constraint called _unbiasedness_[^unbiasedness], which says that $pw(\theta) \geq \alpha$ for all $\theta \in \Theta_a$, and then one can often obtain the UMP unbiased test. This topic is left for more advanced courses.

[^neyman-pearson-level]: Here for simplicity, we assume any choice of $\alpha$ can be exactly achieved, which is the case if the distributions are continuous. The theorem also holds for the discrete case with some technical modifications.
[^hidden-constraint]: Here we have a hidden constraint, which is we need to ensure the level $\alpha$ is correct, i.e. the type I error cannot exceed $\alpha$.
[^unbiasedness]: This is different from the meaning in unbiased estimators.
