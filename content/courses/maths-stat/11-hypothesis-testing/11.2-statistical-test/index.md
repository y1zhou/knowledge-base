---
title: "Statistical Test"
slug: "mathematical-statistics-statistical-test"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics

summary: "Here we introduce the elements of a statistical test, namely null and alternative hypotheses, test statistic, rejection region, and type I and type II errors. We then proceed to large-sample Z-tests and some small-sample tests derived from the small sample CIs."
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
4. More generally, a hypothesis test can be formulated even without introducing a test statistic. One essentially only needs to specify the rejection region on the space of the sample - there is a region $R \subset \mathbb{R}^n$ such that $H_0$ is rejected whenever $(Y_1, \cdots, Y_n) \in R$. The rejection using test statistic $T = t(Y_1, \cdots, Y_n)$ $\in RR$ is a special case  of this, since one can identify $R$ with the preimage of $RR$ under the mapping $t: \mathbb{R}^n \rightarrow R$.



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

# TODO



## Large-sample Z-tests

Now that we know the concepts in a hypothesis test, the important question is how does one construct a hypothesis test. Here we introduce the `Z-test`, which is closely related to the [Z-score]({{< ref "/courses/maths-stat/10-confidence-intervals/index.md#z-score-large-sample-confidence-intervals" >}}).

Suppose we want to test

$$
H_0: \theta = \theta_0 \quad \text{vs.} \quad H_a: \theta \neq \theta_0
$$

where the hypothesized value of interest, $\theta_0$, is regarded as known. Although estimation cannot resolve the decision directly, the following idea is natural:

1. Estimate $\theta$ with an estimator $\hat\theta$;
2. If $\hat\theta$ is far from $\theta_0$, we reject $H_0$.



Now the question becomes *how far is "far"*? We need to take the randomness of $\hat\theta$ into account. If it has a large variance, then we need to think how far it can deviate from $\theta_0$.

**Proposal:** suppose $\hat\sigma$ is the estimated standard error of $\hat\theta$. Let's consider the signed distance between $\theta_0$ and $\hat\theta$ relative to $\hat\sigma$:

$$
Z = \frac{\hat\theta - \theta_0}{\hat\sigma}
$$

Note that unlike the Z-score for constructing the CI, our $Z$ is a statistic because all the quantities can be computed from the sample.

Our proposed rejection rule is we reject $H_0$ if $|Z| > k$, namely when $\hat\theta$ is $k$-times standard error away from $\theta_0$. To choose an appropriate $k$, recall that if $H_0: \theta = \theta_0$ is true, then often when the sample size is large,

$$
Z \overset{\cdot}{\sim} N(0, 1)
$$

so we may choose $k$ using the normal quantile $z_{\alpha/2}$.

To sum up, the proposed rejection rule is


$$
|Z| = \frac{|\hat\theta - \theta_0|}{\hat\sigma} > z_{\frac{\alpha}{2}} \Leftrightarrow \theta_0 > \hat\theta + \hat\sigma z_{\frac{\alpha}{2}} \text{ or } \theta_0 < \hat\theta - \hat\sigma z_{\frac{\alpha}{2}}
$$


The type I error is given by


$$
P_{\theta_0} \left(|Z| > z_{\frac{\alpha}{2}} \right) = \alpha
$$

What we just described is known as the `two-sided Z-test`. There's also one-sided versions, as shown in the table below.


| $H_0$                                         | $H_a$                  | Rejection rule         |
| --------------------------------------------- | ---------------------- | ---------------------- |
| $\theta = \theta_0$                           | $\theta \neq \theta_0$ | $\|Z\| > z_{\alpha/2}$ |
| $\theta = \theta_0$ or $\theta \leq \theta_0$ | $\theta > \theta_0$    | $Z > z_\alpha$         |
| $\theta = \theta_0$ or $\theta \geq \theta_0$ | $\theta < \theta_0$    | $Z < -z_\alpha$        |


In the one-sided Z-tests, one may either have a simple $H_0: \theta = \theta_0$ or a composite $H_0: \theta \leq \theta_0$. Following the calculations in [the example](#example-with-a-normal-distribution), we can show that the rejection rule provides the same correct level for both the simple and the composite cases.

### Bernoulli distribution example
Suppose $Y_i \overset{i.i.d.}{\sim} Bern(\theta)$. Our hypotheses are

$$
H_0: \theta \leq 0.9 \quad \text{vs.} \quad H_a: \theta > 0.9
$$

and we have $n=100$, $\hat\theta = \bar{Y}_n$ and $\alpha = 0.05$. Suppose the observed $\hat\theta = 0.93$. Recall that

$$
Var(\hat\theta) = \frac{\theta(1-\theta)}{n}
$$

Hence the estimated standard error of $\hat\theta$ is

$$
\hat\sigma = \sqrt{\frac{\hat\theta(1-\hat\theta)}{n}} = 0.0255
$$

Then the test statistic is

$$
Z = \frac{0.93-0.9}{0.0255} = 1.18 < z_{0.05} = 1.64
$$

So $Z \notin RR$ and we fail to reject $H_0$.

## Hypothesis test and confidence intervals
You've probably noticed the resemblance between the Z-test and the Z-score confidence interval. In fact, there is a connection between the two. Long story short, whenever you can come up with the CIs, you'd be able to perform the hypothesis test.

Recall we've just discussed that the rejection rule for a two-sided Z-test is
$$
|Z| = \frac{|\hat\theta - \theta_0|}{\hat\sigma} > z_{\frac{\alpha}{2}} \Leftrightarrow \theta_0 > \hat\theta + \hat\sigma z_{\frac{\alpha}{2}} \text{ or } \theta_0 < \hat\theta - \hat\sigma z_{\frac{\alpha}{2}}
$$

On the other hand, a two-sided $(1-\alpha)$-CI based on the Z-score is given by

$$
I = \left[ \hat\theta - \hat\sigma z_{\frac{\alpha}{2}}, \hat\theta + \hat\sigma z_{\frac{\alpha}{2}} \right] = \hat\theta \pm \hat\sigma z_{\frac{\alpha}{2}}
$$

Note how the quantities in the Z-test is exactly the boundaries of the Z-score CI. So the rejection rule can be alternatively expressed as **we reject $H_0$ if $\theta_0 \notin I$**. This means if a CI fails to cover $\theta_0$, then $\hat\theta = \theta_0$ is unlikely.

From this, we may propose a general principle. To `convert a CI to a hypothesis test`, suppose $I$ is a $(1-\alpha)$-CI of $\theta$. To test $H_0: \theta = \theta_0$ vs. $H_a: \theta \neq \theta_0$ at level $\alpha$, we reject $H_0$ if $\theta_0 \notin I$.

The level is correct because under $H_0$, we have

$$
P_{\theta_0}(\theta_0 \notin I) = \alpha
$$

using the definition of [coverage probability]({{< ref "/courses/maths-stat/10-confidence-intervals/index.md#basic-concepts" >}}) of a CI.

### Uniform distribution example
For $Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$, we've derived a one-sided $(1-\alpha)$-CI of the form

$$
I = \left[ \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}} , \infty \right)
$$

where $\hat\theta = \max(Y_1, \cdots, Y_n)$. By the discussion above, we may propose a hypothesis test with the rejection rule as

$$
\theta_0 \notin I \Leftrightarrow \theta_0 < \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}} \Leftrightarrow \hat\theta > \theta_0 (1-\alpha)^\frac{1}{n}
$$

So $H_0: \theta = \theta_0$ or $\theta \leq \theta_0$, and $H_a: \theta > \theta_0$.

It also works the other way around. To `convert a hypothesis test to a CI` (taking the Z-test as an example), the test statistic

$$
Z(\theta_0) = \frac{\hat\theta - \theta_0}{\hat\sigma}
$$

and the rejection region is $\\{z: |z| > z_{\frac{\alpha}{2}}\\}$. Hence,

$$
I = \\{\theta_0 \in \mathbb{R}: Z(\theta_0) \notin RR\\} = \left\\{ \theta_0 \in \mathbb{R}: \left| \frac{\hat\theta - \theta_0}{\hat\sigma} \right| \leq z_{\frac{\alpha}{2}} \right\\} = \left[ \hat\theta - \hat\sigma_{\frac{\alpha}{2}}, \hat\theta + \hat\sigma_{\frac{\alpha}{2}} \right]
$$

which is the Z-score CI. {{<hl>}}The correspondence between CI and HT is one-to-one, although we use $CI \rightarrow HT$ more often.{{</hl>}}

## Some small-sample tests
Recall that we have discussed small sample CIs for means and variances. Based on the correspondence between CI and HT, we can easily develop hypothesis tests based on those CIs.

### One-sample mean
The **assumption** is $Y_i \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$. The hypotheses are

$$
H_0: \mu = \mu_0 \quad vs. \quad H_a: \begin{cases}
  \mu > \mu_0; \\\\
  \mu \neq \mu_0; \\\\
  \mu < \mu_0 \\\\
\end{cases}
$$

The **test statistic**, motivated from [the pivotal quantity for CI]({{<ref "/courses/maths-stat/10-confidence-intervals/index.md#one-sample-mean" >}}), is

$$
T = \frac{\bar{Y}_n - \mu_0}{\hat\sigma / \sqrt{n}} \sim t(n-1) \text{ under }H_0
$$

where $\hat\sigma^2$ is the unbiased estimator of $\sigma^2$:

$$
\hat\sigma^2 = \frac{1}{n-1} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2
$$

The **rejection region** can be found by inverting the CIs:

$$
 RR = \begin{cases}
   \\{t: t > t_\alpha (n-1)\\}; \\\\
   \\{t: |t| > t_\frac{\alpha}{2}(n-1)\\}; \\\\
   \\{t: t < -t_\alpha (n-1)\\}.
 \end{cases}
$$

==Exercise 38==

If we look at $H_0$, is this a simple null or a composite null hypothesis? One might say this is a simple hypothesis since only $\mu_0$ is involved, but in fact this is composite because the unknown $\sigma^2$ is an important part of our model. $H_0$ is essentially

$$
(\mu, \sigma^2) \in \Theta_0 = \\{(\mu_0, u): u > 0\\}
$$

**Def:** an unknown parameter which is not of interest in the hypothesis test is called a `nuisance parameter`, e.g. the $\sigma^2$ above. It's often desirable to consider test statistic $T$ whose distribution does not depend on a nuisance parameter (very similar to the pivotal quantity idea).

The test statistic $T$ above follows a `non-central t-distribution` when $\mu_0$ is **not** the true parameter. This is of interest when analyzing the Type II error/power of the tests, or the type I error under $H_0$: $\mu < \mu_0$ or $\mu > \mu_0$.

### Two-sample means
Our **assumption** is $X_1, \cdots, X_{n_1} \overset{i.i.d.}{\sim} N(\mu_1, \sigma^2)$ and $Y_1, \cdots, Y_{n_2} \overset{i.i.d.}{\sim} N(\mu_2, \sigma^2)$. $X_i$'s are independent of $Y_i$'s. The hypotheses are

$$
H_0: \mu_1 - \mu_2 = \delta_0 \quad vs. \quad H_a: \begin{cases}
  \mu_1 - \mu_2 > \delta_0; \\\\
  \mu_1 - \mu_2 \neq \delta_0; \\\\
  \mu_1 - \mu_2 < \delta_0.
\end{cases}
$$

The **test statistic** again comes from the pivotal quantity in CI:

$$
T = \frac{\bar{X} - \bar{Y} - \delta_0}{\hat\sigma / \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}} \sim t(n_1 + n_2 - 2) \text{ under } H_0
$$

where

$$
\hat\sigma^2 = \frac{\sum_{i=1}^{n_1}(X_i - \bar{X})^2 + \sum_{i=1}^{n_2}(Y_i - \bar{Y})^2}{n_1 + n_2 - 2}
$$

is an unbiased estimator of $\sigma^2$. The $n_1 + n_2 - 2$ is in the denominator because we lost two degrees of freedom when estimating the means. The **rejection regions** are

$$
RR = \begin{cases}
  \\{t: t > t_\alpha(\nu)\\}; \\\\
  \\{t: |t| > t_\frac{\alpha}{2}(\nu)\\}; \\\\
  \\{t: t < -t_\alpha(\nu)\\}.
\end{cases} \quad \text{where } \nu = n_1 + n_2 - 2
$$

When $n_1$ and $n_2$ are very large, the t-distribution becomes really close to the normal distribution, and the Z-score should be used.

### One-sample variance
Our **assumption** is $Y_1, \cdots, Y_n \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$. The hypotheses are

$$
H_0: \sigma^2 = \sigma_0^2 \quad vs. \quad H_a: \begin{cases}
  \sigma^2 > \sigma_0^2; \\\\
  \sigma^2 \neq \sigma_0^2; \\\\
  \sigma^2 < \sigma_0^2.
\end{cases}
$$

The **test statistic** we use is

$$
T = \frac{\sum_{i=1}^n (Y_i - \bar{Y}_n)^2}{\sigma_0^2}
$$

which follows $\chi^2(n-1)$ under $H_0$. The **rejection regions** are

$$
RR = \begin{cases}
  \\{t: t > \chi^2_\alpha (n-1)\\}; \\\\
  \\{t: t > \chi^2_\frac{\alpha}{2} (n-1) \text{ or } t < \chi^2_{1 - \frac{\alpha}{2}}(n-1)\\}; \\\\
  \\{t: t < \chi^2_{1 - \alpha}(n-1)\\}.
\end{cases}
$$

where $\chi^2_\alpha (\nu)$ denotes the $(1-\alpha)$-quantile of $\chi^2(\nu)$.

> There's also a two-sample variance test involving the F-distribution.
