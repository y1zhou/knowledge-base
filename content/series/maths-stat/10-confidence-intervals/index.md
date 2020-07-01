---
title: "Confidence Intervals"
slug: "mathematical-statistics-confidence-intervals"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - CI

summary: "Confidence intervals and methods of contructing them."
date: 2020-02-08T09:57:16-04:00
toc: true
weight: 100
---

## Basic concepts

When we talk about estimation, everything that comes in should be known, i.e. computable from the sample. This idea is also reflected in the definition of confidence intervals.

**Def:** A `confidence interval` (CI) is given by one of the following three forms:
$$
I = \begin{cases}
	[L, U], & \text{two-sided}, \\\\
	(-\infty, U], \\\\
	[L, \infty), & \text{one-sided}
\end{cases}
$$
where the lower-bound $L$ and upper-bound $U$ are *statistics*. **Never** put unknown quantities in CIs!

Most often, $L$ and $U$ have continuous distributions, so taking closed or open intervals doesn't matter. Another very important remark is since the sample is (in this course at least) random, the CIs are also random in general. On the other hand, when a CI is reported from the data, it's *not* random anymore.

For example, suppose $Y \sim N(\theta, 1)$ where $\theta \in \mathbb{R}$.

- $[Y - 2, Y + 2]$ is a CI because both $L$ and $U$ are known.
- $[-2, 2]$ is a CI.
- $[\theta - 2, \infty)$ is not a confidence interval because $\theta$ is unknown.
- $[Y-2, \infty)$ is a CI.

**Def:** Assume $\theta$ is the unknown true parameter, and $I$ is the CI, either one-sided or two-sided. $\alpha = P(\theta \notin I) = 1 - P(\theta \in I)$ is called the `confidence level` or confidence coefficient. $P(\theta \in I) = 1 - \alpha$ is called the `coverage probability` or the `confidence coefficient`. We also say $I$ is a $(1-\alpha)$-CI. Often $\alpha$ is chosen as $0.05$, which gives a 95% confidence interval.

Note here $I$ is random, and $\theta$ is not! The randomness is in the interval, not in the parameter, at least from the frequentist's point of view. In Bayesian statistics, $\theta$ is also random and they talk more about credible intervals.

**Remark:** Ideally, given the coverage probability, we want the CI to be as small as possible. We will not address the *optimality* of CI in this course. Later we'll introduce the optimality for hypothesis testing, which in a way translates to the optimality of CI.

### Normal distribution example

Suppose we have $Y \sim N(\theta, 1)$. Propose a one-sided CI for $\theta$.

If we have a CI of $[Y, \infty)$, the coverage probability is 0.5 because it's equally likely $Y$ is smaller or larger than the true $\theta$.

**Def:** for $x \in (0, 1)$, the `tail quantile` $z_x$ of a standard normal  $Z \sim N(0, 1)$ is defined as the number satisfying
$$
P(Z > z_x) = x = P(Z < -z_x) = P(Z < z_{1-x})
$$
We claim that $[Y-z_\alpha, \infty)$ is a $(1-\alpha)$-CI for $\theta$. To check this, first of all note that $Y-Z_\alpha$ is still a statistic.
$$
P(\theta \in [Y-z_\alpha, \infty)) = P(Y - \theta \leq z_\alpha) = P(Z \leq z_\alpha) = 1-\alpha
$$
Similarly, one can check that $(-\infty, Y+z_\alpha]$ is also a $(1-\alpha)$-CI.

### Uniform distribution example

Suppose $Y_1, \cdots, Y_n$ are i.i.d. from $Unif(0, \theta)$. We want to find a one-sided CI of form $[L, \infty)$.

We know that $\hat\theta = \max(Y_1, \cdots, Y_n)$ is the MLE. We claim that
$$
\left[ \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}, \infty \right)
$$
is the $(1-\alpha)$-CI for $\theta$. $\frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}$ is obviously a statistic, and now we want to check its coverage probability.
$$
\begin{gather*}
	P\left(\theta \geq \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}\right) = P\left( \hat\theta \leq \theta(1-\alpha)^\frac{1}{n} \right) = P\left( Y_1 \leq \theta(1-\alpha)^\frac{1}{n} \right)^n \\\\
	F(y) = \frac{y}{\theta}, \quad 0 \leq y \leq 1 \\\\
	P\left(\theta \geq \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}\right) = \left(\frac{\theta(1-\alpha)^\frac{1}{n}}{\theta}\right)^n = 1-\alpha
\end{gather*}
$$
If $\alpha = 0.05$, suppose from the data we got $\hat\theta = 0.9$ and the sample size $n = 5$. Plugging these values in and we can find the observed CI to be $[0.9234, \infty)$.

**Warning:** If a CI $I_0$ is observed from data, it **never** makes sense to write $P(\theta \in I_0) = 1-\alpha$. For example,
$$
P(\theta \in [0.9234, \infty)) = 0.95
$$
because $\theta$ is unknown but *fixed*. The "probability" should be either 0 or 1.



### Two-sided cases

Suppose $L$ and $U$ are statistics based on some sample. We know that

- $L < U$,
- $(-\infty, U]$ is a $(1-\alpha_1)$-CI for $\theta$, and
- $[L, \infty)$ is a $(1-\alpha_2)$-CI for $\theta$,

then $[L, U]$ is a $(1-\alpha_1 - \alpha_2)$-CI for $\theta$. To prove this, define events $A = \\{\theta \leq U\\}$ and $B = \\{\theta \geq L\\}$. We know $P(A) = 1-\alpha_1$ and $P(B) = 1 - \alpha_2$.
$$
\begin{aligned}
	P(L \leq \theta \leq U) &= P(A \cap B) = P(A) + P(B) - P(A \cup B) \\\\
	&= (1-\alpha_1) + (1-\alpha_2) - 1 \\\\
	&= 1 - \alpha_1 - \alpha_2
\end{aligned}
$$
using the inclusion-exclusion principle, and $P(A \cup B) = 1$ because $L < U$.

If $Y \sim N(\theta, 1)$, we know that $[Y - z_{\frac{\alpha}{2}}, \infty)$ and $(-\infty, Y + z_{\frac{\alpha}{2}}]$ are both $(1-\frac{\alpha}{2})$-CIs. By the proposition above, $\left[Y - z_{\frac{\alpha}{2}}, Y -+z_{\frac{\alpha}{2}} \right]$ is a $(1-\alpha)$-CI.

## Construction of confidence intervals

A quantity $T = f(\theta, Y_1, \cdots, Y_n)$ is said to be `pivotal` if, assuming $\theta$ is true, the *distribution* of $T$ does not depend on $\theta$.

**Remark:** A pivotal quantity is not a statistic in general. Although its functional form depends on $\theta$, its distribution doesn't. For example, $Y \sim N(\theta, 1)$ where $\theta$ is unknown. $T = Y - \theta$ is pivotal because the distribution is $N(0, 1)$, but it's not a statistic since it's not computable from the sample.

### Uniform distribution example

$Y_1, \cdots, Y_n \sim Unif(0, \theta)$. $\hat\theta = \max(Y_1, \cdots, Y_n)$. We claim that $T = \frac{\hat\theta}{\theta}$ is pivotal.
$$
T = \max\left(\frac{Y_1}{\theta}, \cdots, \frac{Y_n}{\theta}\right)
$$
We can observe that $\frac{Y_i}{\theta}$ are i.i.d. from $Unif(0, 1)$. This can be checked using the CDF.
$$
P(T \leq y) = P\left(\frac{Y_1}{\theta} \leq y\right)^n = y^n, \quad 0 \leq y \leq 1
$$
which doesn't depend on $\theta$.

### Pivotal method

We now introduce the pivotal method for constructing one-sided $(1-\alpha)$-CIs. The uniform distribution above is used again as an example for explaining the procedure.

1. Construct a pivotal quantity involving $\theta$, the parameter of interest. $T = \frac{\hat\theta}{\theta}$.
2. Setup an equation of either of the following forms:
   1. $P(T \leq y) = 1 - \alpha$. We got $P(T \leq y) = y^n = 1-\alpha$, and solving for $y$ yields $y = (1-\alpha)^\frac{1}{n}$.
   2. $P(T \geq y) = 1 - \alpha$.
3. Solve the inequality $T \leq y$ or $T \geq y$ for $\theta$. $T = \frac{\hat\theta}{\theta} = (1-\alpha)^\frac{1}{n}$. Solving for $\theta$ yields $\theta \geq \frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}$.
4. Conclude that $\left[\frac{\hat\theta}{(1-\alpha)^\frac{1}{n}}, \infty \right)$ is a $(1-\alpha)$-CI for $\theta$. Note that we got the lower bound because we chose $P(T \leq y)$ in step 2.

A natural question that arises is why does $T$ need to be pivotal. The answer is if $T$ is not pivotal, the solution for $y$ in step 2 will depend on $\theta$, which means the CI will depend on $\theta$.



### Normal distribution example

Suppose $Y \sim N(\theta, 1)$, and we want to find the two-sided CI for $\theta$.

We know that $Y - \theta \sim N(0, 1)$ is pivotal. 
$$
\begin{gather*}
	P(T \leq y_1) = 1 - \frac{\alpha}{2} && P(T \geq y_2) = 1 - \frac{\alpha}{2} \\\\
	y_1 = z_{\frac{\alpha}{2}} && y_2 = z_{1 - \frac{\alpha}{2}} = -z_{\frac{\alpha}{2}} \\\\
	T \leq y_1 \Leftrightarrow Y - \theta \leq z_{\frac{\alpha}{2}} \Leftrightarrow \theta \geq Y - z_{\frac{\alpha}{2}} && T \geq y_2 \Leftrightarrow \theta \leq Y + z_{\frac{\alpha}{2}}
\end{gather*}
$$
From this we can see that both $\left[Y - z_{\frac{\alpha}{2}}, \infty \right)$ and $\left(-\infty, Y + z_{\frac{\alpha}{2}} \right]$ are $(1-\frac{\alpha}{2})$-CIs, so $\left[Y - z_{\frac{\alpha}{2}}, Y + z_{\frac{\alpha}{2}} \right]$ is a $(1-\alpha)$-CI for $\theta$.

**Remark:** As in the example above, we often (by default) choose $y_1$ and $y_2$ such that $P(T \leq y_1) = P(T \geq y_2)$. This is however not the only option for constructing CIs. In fact, if we dive in the optimality of CIs there are multiple ways of doing it, such as finding the narrowest CI.



## Z-score (large-sample) confidence intervals

The most challenging step in the pivotal method is to find the pivotal quantity (step 1). Fortunately, for many common situations when the sample size is large, one type of pivotal quantity called the `Z-score` can be used.

### Motivating discussion

Suppose $Y_i$ are i.i.d. samples from some distribution with mean $E[Y_i] = \mu$ and variance $Var(Y_i) = s^2$, where $\mu$ is unknown and $s^2$ is typically unknown. Our goal is to construct a $(1-\alpha)$-CI for $\mu$.

Recall that $E \left[\bar{Y}_n \right] = \mu$, and $Var\left(\bar{Y}_n \right) = \frac{s^2}{n}$. If we try to standardize $\bar{Y}_n$,


$$
\frac{\bar{Y}_n - E[\bar{Y}_n]}{\sqrt{Var(\bar{Y}_n)}} = \frac{\bar{Y}_n - \mu}{s / \sqrt{n}} \overset{d}{\approx} N(0, 1)
$$


holds by the CLT. Note that if $s^2$ is unknown, we can replace it with $S^2$, which is a consistent estimator of $s^2$. For example, $S^2 = \frac{1}{n-1}\sum_{i=1}^n (Y_i - \bar{Y}_n)^2$. Then,
$$
Z = \frac{\bar{Y}_n - \mu}{S / \sqrt{n}} \overset{d}{\approx} N(0, 1)
$$

> $\overset{d}{\approx}$ means approximately equal in distribution.

We can see that $Z$ is (approximately) pivotal. Note that $P(-z_{\frac{\alpha}{2}} \leq Z \leq z_{\frac{\alpha}{2}}) \approx 1 - \alpha$. If we solve this for $\mu$, we can get
$$
\bar{Y}_n - z_{\frac{\alpha}{2}} \cdot \frac{S}{\sqrt{n}} \leq \mu \leq \bar{Y}_n + z_{\frac{\alpha}{2}} \cdot \frac{S}{\sqrt{n}}
$$
which is the two-sided CI $\left[ \bar{Y}_n - z_{\frac{\alpha}{2}} \cdot \frac{S}{\sqrt{n}}, \bar{Y}_n + z_{\frac{\alpha}{2}} \cdot \frac{S}{\sqrt{n}} \right]$.

### Definition

In fact, this idea can be generalized. Suppose $\hat\theta$ is an estimator of $\theta$, and $\hat\sigma^2$ is a consistent estimator of $\sigma^2$, the variance of $\hat\theta$. The `Z-score` is defined as
$$
Z = \frac{\hat\theta - \theta}{\hat\sigma}
$$
and $\hat\theta$ is asymptotically unbiased. Note that we can plug $\sigma$ in directly if it's known. Although $\hat\theta$ may not be a sample mean, it turns out often $Z \overset{d}{\approx} N(0, 1)$ when the sample size is large.

For example, if $\hat\theta$ is the MLE of $\theta$, then typically $\hat\theta \overset{d}{\approx} N(\theta, \frac{1}{nI(\theta)})$ where $I(\theta)$ is the Fisher information. Based on this fact, we may propose
$$
Z = \frac{\hat\theta - \theta}{\sqrt{nI(\hat\theta)}}
$$
as the Z-score. We shall assume that a Z-score follows $N(0, 1)$ omitting justification in this course.

Using the Z-score as a pivotal quantity, a $(1-\alpha)$-CI can be constructed as $\left[ \hat\theta - \hat\sigma z_{\frac{\alpha}{2}}, \hat\theta + \hat\sigma z_{\frac{\alpha}{2}} \right]$.

### Example

Suppose we have the data for two players: 

| Player | Shooting % | Number of shots |
| ------ | ---------- | --------------- |
| 1      | 40%        | 200             |
| 2      | 50%        | 100             |

Our goal is to construct a 95% CI for $p_1$, $p_2$ and $\theta = p_1 - p_2$.

Suppose $X_i$ represents the outcome for Player 1. We have $X_i \overset{i.i.d.}{\sim} Bern(p_1)$ and similarly $Y_i \overset{i.i.d.}{\sim} Bern(p_2)$. We assume that $\\{X_i\\}$ and $\\{Y_i\\}$ are independent. Recall that the Z-score CI is given as $\hat\theta \pm z_{\frac{\alpha}{2}} \hat\sigma$.


$$
\begin{aligned}
	\hat{p}_1 &= \bar{X}\_{n_1} \xlongequal{\text{obs. as}} 0.4 \\\\
	\hat{p}_2 &= \bar{Y}\_{n_2} = 0.5 \\\\
	\hat\theta &= \hat{p}_1 - \hat{p}_2 = -0.1
\end{aligned}
$$


Now we need to find an estimation for the standard deviation.
$$
\begin{gather*}
	s_1^2 = Var(X_1) = p_1(1-p_1) \\\\
	s_2^2 = Var(Y_1) = p_2(1-p_2) \\\\
	Var(\hat{p}_1) = \frac{s_1^2}{n_1} \\\\
	Var(\hat{p}_2) = \frac{s_2^2}{n_2} \\\\
	Var(\hat\theta) = \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2}
\end{gather*}
$$
Note that $p_1$ and $p_2$ are unknown, but we can plug in their estimates to get
$$
\begin{gather*}
	\hat\sigma_1 = \sqrt{\frac{\hat{p}_1(1 - \hat{p}_1)}{n_1}} = \sqrt{\frac{0.4 \times 0.6}{200}} \approx 0.03464 \\\\
	\hat\sigma_2 = \sqrt{\frac{\hat{p}_2(1 - \hat{p}_2)}{n_2}} = \sqrt{\frac{0.5 \times 0.5}{100}} = 0.05 \\\\
	\hat\sigma = \sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}} \approx 0.06083
\end{gather*}
$$
Since $z_{0.025} = 1.96$, the 95%-CI for $p_1$, $p_2$ and $\theta$ are respectively given by
$$
\begin{gather*}
	0.4 \pm 0.03464 \times 1.96 = [0.33, 0.47] \\\\
	0.5 \pm 0.05 \times 1.96 = [0.40, 0.60] \\\\
	-0.1 \pm 0.06083 \times 1.96 = [-0.22, 0.02]
\end{gather*}
$$
The CI for $\theta$ covers $0$, so we can't say with confidence that Player 2 is better. We can also see that the CI for $p_1$ is the shortest, and this is because of different sample sizes since $\hat\sigma \approx s.e.(\hat\theta)$.

The quantity $z_\frac{\alpha}{2} \hat\sigma$ is known as the `margin of error` because it's half the length of the Z-score CI. This quantity can be used for selecting sample size. See Section 8.7 in the textbook for a more detailed explanation.



## Some small-sample CIs

When the sample size is small, one cannot resort to the Z-score CI. To compensate the lack of information due to the small sample size, we shall make a strong assumption: the sample follows a normal distribution. This normal assumption is supported by the CLT:
$$
X = Y_1 + \cdots + Y_n \overset{d}{\approx} N(E[X], Var(X))
$$


> This assumption is only valid when the distribution is continuous.

Assume $Z, Z_1, Z_2, \cdots, Z_n \overset{i.i.d.}{\sim}N(0, 1)$.

A $\chi^2$-distribution with $\nu$ degrees of freedom is the distribution of $W = Z_1^2 + \cdots + Z_\nu^2$.

The $t$-distribution with $\nu$ d.f. is the distribution of $T = \frac{Z}{\sqrt{W/\nu}}$.

- $t(\nu)$ has a similar shape to $N(0, 1)$.
- When $\nu$ is small ($\nu \leq 30$), the $t(\nu)$ distribution has noticably heavier tails than $N(0, 1)$.
- As $\nu \rightarrow \infty$, $t(\nu)$ becomes $N(0, 1)$.

### One-sample mean

Assume $Y_i \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$, and both parameters are unknown. Our goal is to construct a CI for $\mu$.

We know that


$$
T_n = \frac{\bar{Y}_n - \mu}{\hat\sigma_n / \sqrt{n}} \sim t(n-1)
$$


where $\hat\sigma_n = \frac{1}{n-1} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2$. $T_n$ is a pivotal quantity since its distribution is free of $\mu$ or $\sigma^2$. Using the pivotal method, a two-sided $(1-\alpha)$-CI for $\mu$ is given as
$$
\bar{Y}_n \pm t_{\frac{\alpha}{2}}(n-1) \frac{\hat\sigma_n}{\sqrt{n}}
$$
where $t_{\frac{\alpha}{2}}(n-1)$ is the $\alpha/2$ tail quantile of a $t$-distribution with $n-1$ degrees of freedom.



### Two-sample mean

Assume $X_i \overset{i.i.d.}{\sim} N(\mu_1, \sigma^2)$ and $Y_i \overset{i.i.d.}{\sim} N(\mu_2, \sigma^2)$. The sample sizes are $n_1$ and $n_2$, respectively, and $X_i$ and $Y_i$ are independent. Our goal is to construct a CI for $\mu_1 - \mu_2$.

Recall that


$$
s.e.(\bar{X}\_{n_1} - \bar{X}\_{n_2}) = \sqrt{\frac{\sigma^2}{n_1} + \frac{\sigma^2}{n_2}} = \sigma\sqrt{\frac{1}{n_1} + \frac{1}{n_2}}
$$


The pivotal quantity is given by


$$
T = \frac{(\bar{X}\_{n_1} - \bar{X}\_{n_2}) - (\mu_1 - \mu_2)}{\hat\sigma_{X, Y} \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}}
$$


where


$$
\hat\sigma_{X, Y}^2 = \frac{\sum\_{i=1}^{n_1} (X_i - \bar{X}_{n_1})^2 + \sum\_{i=1}^{n_2} (Y_i - \bar{Y}_{n_2})^2}{n_1 + n_2 - 2}
$$




The $(n_1 + n_2 - 2)$ makes $\hat\sigma_{X, Y}$ unbiased. $T \sim t(n_1 + n_2 - 2)$ so $T$ is pivotal. Applying the pivotal method yields the CI
$$
(\bar{X}_{n_1} - \bar{Y}_{n_2}) \pm t_{\frac{\alpha}{2}}(n_1 + n_2 - 2) \hat\sigma_{X, Y}\sqrt{\frac{1}{n_1} + \frac{1}{n_2}}
$$


If the variance of the two populations are different, refer to `Welch's t-statistic`.

### CI for $\sigma^2$

We assume $Y_i \overset{i.i.d.}{\sim} N(\mu, \sigma^2)$, and our goal is to construct a CI for $\sigma^2$.

**Theorem:** Let $Y_1, \cdots, Y_n$ be a random sample from a normal distribution with mean $\mu$ and variance $\sigma^2$. Then


$$
\frac{(n-1)S^2}{\sigma^2} = \frac{1}{\sigma^2}\sum_{i=1}^n (Y_i - \bar{Y})^2
$$

has a $\chi^2$ distribution with $(n-1)$ d.f. Also, $\bar{Y}$ and $S^2$ are independent random variables. By this theorem,

$$
T = \frac{\sum_{i=1}^n (Y_i - \bar{Y}_n)^2}{\sigma^2} = \frac{(n-1)\hat\sigma_n^2}{\sigma^2} \sim \chi^2(n-1)
$$

and $T$ is once again a pivotal quantity. We choose $x < y$ such that $P(T \leq x) = P(T \geq y) = \frac{\alpha}{2}$. We're doing this because the $\chi^2$ distribution is asymmetric. Denote such $x, y$ as

$$
\chi^2_{(1 - \frac{\alpha}{2})}(n-1) \text{  and  } \chi^2_{\frac{\alpha}{2}}(n-1)
$$

A $(1-\alpha)$-CI is given by

$$
\left[ \frac{(n-1)\hat\sigma_n^2}{\chi^2_{\alpha/2}(n-1)}, \frac{(n-1)\hat\sigma_n^2}{\chi^2_{(1-\alpha/2)}(n-1)} \right]
$$

**Remark:** the validity of the CIs introduced in this section depends on the normality assumption. Empirically, people find that these CIs are robust against moderate departure from the normality assumption.
