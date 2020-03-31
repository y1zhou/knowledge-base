---
title: "The Method of Moments"
slug: "mathematical-statistics-method-of-moments"

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Estimation

summary: "A fairly simple method of constructing estimators that's not often used now."
date: 2020-01-28T10:46:18-04:00
toc: false
type: docs  # Do not modify.
weight: 130

menu:
  maths-stat:
    name: Method of Moments
    parent: Estimation
    weight: 30
---

**Def:** Suppose $Y_1, \cdots, Y_n$ are i.i.d. samples. The $k$-th `population moment` is
$$
\mu_k = E\left[ Y_i^k\right ].
$$

The $k$-th `sample moment` is defined as

$$
m_k = \frac{1}{n}\sum_{i=1}^n Y_i^k, k = 1, 2, 3, \cdots
$$

**Remark:** All moments $m_1, m_2, \cdots$ are statistics. If $E\left[ |Y_i|^k\right ] < \infty$, then by the law of large numbers,

$$
m_k = m_{k, n} = \frac{1}{n}\sum_{i=1}^n Y_i^k \xrightarrow{P} E\left[ Y_i^k\right ] = \mu_k
$$

If our goal is to find the population moments, this is good enough. For example, $Y_i \overset{i.i.d.}{\sim} Unif(0, \theta)$. Recall that $2\bar{Y}_n$ was proposed as an estimator of $\theta$. This time we use the method of moments to derive the estimator.

$$
\begin{aligned}
	\mu_1 &= E[Y_i] = \frac{\theta}{2} \\
	m_1 &= \frac{1}{n}\sum_{i=1}^n Y_i = \bar{Y}_n
\end{aligned}
$$

Letting $\mu_1 = m_1$, we get

$$
\frac{\theta}{2} = \bar{Y}_n \Rightarrow \hat\theta = 2\bar{Y}_n
$$

## Procedure

$Y_i$ are i.i.d. samples from a distribution with $r$ unknown parameters $\theta_1, \cdots, \theta_r$.

1. Compute population moments.

$$
\begin{cases}
	\mu_1 = \mu_1(\theta_1, \theta_2, \cdots, \theta_r) = E\left[ Y_1\right ] \\\\
	\mu_2 = E\left[ Y_1^2\right ] \\\\
	\vdots \\\\
	\mu_r = E\left[ Y_1^r\right ]
\end{cases}
$$

2. Write down a system of equations in the following way:

$$
\begin{equation} \label{eq:sample-moments}
\begin{cases}
	\mu_1(\theta_1, \cdots, \theta_r) = m_1 \\\\
	\mu_2(\theta_1, \cdots, \theta_r) = m_2 \\\\
	\vdots \\\\
	\mu_r(\theta_1, \cdots, \theta_r) = m_r
\end{cases}
\end{equation}
$$

3. Solve $\eqref{eq:sample-moments}$ to get the estimators of $\theta_1, \cdots, \theta_r$. That is, express $\theta_1, \cdots, \theta_r$ in terms of statistics $m_1, \cdots, m_r$.

**Remark:** The method of moments often yield *consistent* estimators. The system of equations $\eqref{eq:sample-moments}$ can be rewritten in vector form as
$$
\begin{gather*}
	\vec{\theta} = (\theta_1, \cdots, \theta_r), \qquad \vec{m} = (m_1, \cdots, m_r) \\\\
	\vec{\mu}(\vec{\theta}) = \bigg( \mu_1(\theta_1, \cdots, \theta_r), \cdots, \mu_r(\theta_1, \cdots, \theta_r) \bigg) \\\\
	\eqref{eq:sample-moments} \Leftrightarrow \vec{\mu}(\vec{\theta}) = \vec{m}
\end{gather*}
$$

So what does solving the system of equations actually mean? Suppose there exists an inverse of $\vec{\mu}$, denoted as $\vec{\mu}^{-1}$ continuous. Recall that by the law of large numbers, $m_k \xrightarrow{P} \mu_k$ as $n \rightarrow \infty$.

$$
\hat\theta = \vec{\mu}^{-1}(\vec{m}) = \vec{\mu}^{-1}((m_1, \cdots, m_r)) \xrightarrow{P} \vec{\mu}^{-1}(\mu_1, \cdots, \mu_r) = (\theta_1, \cdots, \theta_r)
$$

As an example, suppose $Y_i \overset{i.i.d.}{\sim} Gamma(\alpha, \beta)$. We know that $E[Y_i] = \alpha\beta$ and $Var(Y_i) = \alpha\beta_2$. To find estimators for $\alpha$ and $\beta$, the population moments are

$$
\begin{gather*}
	\mu_1(\alpha, \beta) = \alpha\beta \\\\
	\mu_2(\alpha, \beta) = Var(Y_i) + E[Y_i]^2 = \alpha\beta_2 + \alpha^2\beta_2
\end{gather*}
$$

Set the system of equations as

$$
\begin{cases}
	m_1 = \alpha\beta \\\\
	m_2 = \alpha^2\beta^2 + \alpha\beta_2
\end{cases}
$$

Substitute the first equation into the second one and we get

$$
m_2 = m_1^2 + m_1\beta \Rightarrow \beta = \frac{m_2 - m_1^2}{m_1}
$$

Putting this back into the first equation

$$
\alpha = \frac{m_1}{\beta} = \frac{m_1^2}{m_2 - m_1^2}
$$

And we can get our estimates:
$$
\begin{aligned}
	\hat\alpha &= \frac{m_1^2}{m_2 - m_1^2} \xrightarrow{P} \frac{\mu_1^2}{\mu_2 - \mu_1^2} = \alpha \\\\
	\hat\beta &= \frac{m_2 - m_1^2}{m_1} \xrightarrow{P} \frac{\mu_2 - \mu_1^2}{\mu_1} = \beta
\end{aligned}
$$

This way of constructing estimators is [fairly simple](https://en.wikipedia.org/wiki/Method_of_moments_(statistics)), but is not used much nowadays because we have better methods, as introduced in the [next chapter]({{< ref "/courses/maths-stat/9-estimation-under-parametric-models/9.1-maximum-likelihood-estimator/index.md" >}}).
