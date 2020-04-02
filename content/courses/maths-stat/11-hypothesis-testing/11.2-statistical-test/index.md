---
title: "Statistical Test"
slug: "mathematical-statistics-statistical-test"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics

summary: ""
date: 2020-04-01T16:55:50-04:00
toc: true
type: docs  # Do not modify.
weight: 190

menu:
  maths-stat:
    name: Statistical Test
    parent: Hypothesis Testing
    weight: 20
---

The imbalanced severity of false decision leads to the concept of hypothesis testing. In a statistical decision where mistakes are imbalanced, we first guard against the mistake of more severity, say to control the probability to be smaller than 0.05, and then do our best to minimize the probability of the other mistake.

## Elements of a statistical test

**Def:** `Statistical hypothesis test` (HT) is a statistical decision process which has the following elements:

- Null hypothesis $H_0: \theta \in \Theta_0$.
- Alternative hypothesis $H_a: \theta \in \Theta_a$ where $\Theta_0$ and $\Theta_a$ are disjoint, and $\Theta_0 \cup \Theta_a = \Theta$. We also have that the severity of falsely accepting $H_0 <$ the severity of falsely accepting $H_a$.
  - This is also how we decide which hypothesis is $H_0$.
- Test statistic $T$.
- Rejection region $RR$.
- Decision rule: if $T \in RR$, we reject $H_0$ or accept $H_a$; otherwise we fail to reject $H_0$.

> When $T \notin RR$, it's often preferable to say "fail to reject $H_0$" than "accept $H_0$", since there's a built-in bias towards guarding $H_0$. In other words, we're being very conservative about $H_0$.



### Type I and type II errors

Suppose we have i.i.d. sample $Y_1, \cdots, Y_n \sim N(\theta, 1)$. Our hypotheses are


$$
H_0: \theta \leq 0 \quad \text{vs.} \quad H_a: \theta > 0
$$


The parameter spaces are $\Theta_0 = (-\infty, 0]$ and $\Theta_a = (0, \infty)$. We propose the sample mean as a test statistic $T = \bar{Y}_n$, which is the estimator of $\theta$.

Intuitively, when $\theta \leq 0$, $\bar{Y}_n$ is unlikely to be large. Because of this, we propose $RR = [c, \infty)$ where $c$ is a to be determined threshold, and we reject $H_0$ if and only if $T \geq c$.

**Def:** In a hypotheses test, the `type I error` is when $H_0$ is true but we choose $H_a$ (falsely reject $H_0$), whose probability is denoted as


$$
\alpha(\theta) = P_\theta(T \in RR), \quad \theta \in \Theta_0
$$
`Type II error` is when $H_a$ is true but we choose $H_0$ (falsely accept $H_0$), and the probability is denoted as


$$
\beta(\theta) = P_\theta(T \notin RR), \quad \theta \in \theta_a
$$


The largest possible type I error probability is called the `level` of the test:


$$
\alpha = \max \Big(\alpha(\theta): \theta \in \Theta_0 \Big)
$$


> Strictly speaking, we should use the supremum instead of the maximum here.



We shall choose the size of $RR$ so that $\alpha$ matches a given level, say 0.05. The probability of accepting $H_a$ when $H_a$ is true is called the `power` of the test. Namely,


$$
pw(\theta) = P_\theta(T \in RR) = 1 - \beta(\theta), \quad \theta \in \Theta_a
$$


#### Remarks

1. In general, both $\alpha(\theta)$ and $\beta(\theta)$ are functions of $\theta$, but defined on different domains $\Theta_0$ and $\Theta_a$, respectively.

2. If $\Theta_0$ consists of a single element $\theta_0$, we say that we have a `simple null hypothesis`. Otherwise, we say we have a `composite null hypothesis`. Similarly we have simple/composite alternative hypotheses. In particular when $\Theta_0 = \\{\theta_0\\}$, the level is directly given as


$$
\alpha = \alpha(\theta_0) = P_{\theta_0}(T \in RR)
$$


3. Power, the probability of correctly accepting $H_a$, is often interpreted as the ability to detect $H_a$. {{<hl>}} $pw(\theta)$ typically increases when the sample size increases, or when $\theta \in \Theta_a$ moves away from the boundary between $\Theta_0$ and $\Theta_a$. {{</hl>}}
4. More generally, a hypothesis test can be formulated even without introducing a test statistic. One essentially only needs to specify the rejection region on the space of the sample - there is a region $R \sub \mathbb{R}^n$ such that $H_0$ is rejected whenever $(Y_1, \cdots, Y_n) \in R$. The rejection using test statistic $T = t(Y_1, \cdots, Y_n)$ $\in RR$ is a special case  of this, since one can identify $R$ with the preimage of $RR$ under the mapping $t: \mathbb{R}^n \rightarrow R$.



### Example with a normal distribution

Following the example where $Y_i \overset{i.i.d.}{\sim} N(\theta, 1)$, $H_0: \theta \leq 0 $, $H_a: \theta > 0$, $T = \bar{Y}_n$ and $RR = [c, \infty)$. Note here we fixed $RR$, but in reality there's some flexibility where we decide the size of $RR$ based on the desired level of $\alpha$.

When $\theta$ is the true parameter, we know that


$$
\bar{Y}_n \sim N\left(\theta, \frac{1}{n} \right) \Leftrightarrow \frac{\bar{Y}_n - \theta}{\sqrt{1/n}} \sim N(0, 1)
$$


#### Type I error

When $\theta \in \Theta_0 = (-\infty, 0]$, the **type I error** is given by


$$
\begin{aligned}
  \alpha(\theta) &= P_\theta(T \in RR) = P_\theta(\bar{Y}_n \geq c) \\\\
  &= P\left( \frac{\bar{Y}_n - \theta}{\sqrt{1/n}} \geq \frac{c - \theta}{\sqrt{1/n}} \right) \\\\
  &= 1 - \Phi \left(\sqrt{n}(c-\theta) \right)
\end{aligned}
$$


where $\Phi$ is the CDF of $N(0, 1)$. The **level** is


$$
\alpha = \max_{\theta \leq 0} \Big( 1 - \Phi \left(\sqrt{n}(c-\theta) \right) \Big)
$$


Note that since the CDF $\Phi$ is monotonically increasing, the function $1 - \Phi \left(\sqrt{n}(c-\theta) \right)$ is also monotonically increasing in $\theta$, so the maximum is achieved at the boundary $\theta = 0$:


$$
\alpha = 1 - \Phi(c\sqrt{n})
$$


So $c\sqrt{n} = z_{\alpha}$, and we need to choose


$$
c = \frac{z_\alpha}{\sqrt{n}}
$$


where $z_\alpha$ is the $(1-\alpha)$-quantile of $N(0, 1)$.



#### Type II error

When $\theta \in \Theta_a = (0, \infty)$,


$$
\begin{aligned}
  \beta(\theta) &= P_\theta(T \notin RR) = P_\theta(\bar{Y}_n < c) \\\\
  &= P_\theta \left( \sqrt{n}(\bar{Y}_n - \theta) < \sqrt{n}(c-\theta) \right) \\\\
  &= \Phi \left(\sqrt{n}(c-\theta) \right)
\end{aligned}
$$



When we plug in the previously worked out $c = z_\alpha / \sqrt{n}$,


$$
\beta(\theta) = \Phi \left( z_\alpha - \sqrt{n}\theta \right), \quad \theta > 0
$$


The **power** is given by


$$
pw(\theta) = 1 - \Phi \left( z_\alpha - \sqrt{n}\theta \right), \quad \theta > 0
$$


Note how the power increases as $\theta$ increases (moves away from the boundary 0), and increases as $n$ increases.



### Two-sided hypothesis

Suppose $Y_i$ are i.i.d. $N(\theta, 1)$ as above. This time we want to test $H_0: \theta = 0$ vs. $H_a: \theta \neq 0$. The test statistic $T$ is still $\bar{Y}_n$, but $RR$ now takes the form


$$
(-\infty, -c] \cup [c, \infty)
$$


where $c$ is a threshold chosen to make the test at level $\alpha$. The calculation of the type I, II error probabilities and the power is given as follows.

