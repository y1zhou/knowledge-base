---
title: "Correlation and Concordance"
slug: "nonparametric-methods-correlation-and-concordance"
categories:
  - Nonmarametric Methods
tags:
  - Nonmarametric Methods
  - Statistics
  - Correlation

summary: "Measures for the strength of relationships between variables (two or more). The Spearman rank correlation coefficient, Kendall's tau and Kendall's W are introduced."
date: 2019-05-05T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 80

menu:
  nonparam-stat:
    name: Correlation and Concordance
    parent: Association Analysis
    weight: 10
---

We've been focusing on location inference for quite a while. There's of course other inferences in the field, and what we often want are measures that summarize the strength of relationships between variables, namely the strength of association or dependence.

Recall that we have the "classical" Pearson correlation coefficient between two random variables $X$ and $Y$. It's a measure of **linear association**. Inference for $\rho$ (population parameter) based on $r$ (sample value) has an assumption of `bivariate normality`, i.e. $X$ and $Y$ are jointly normally distributed.

Can we be more general and relax away the normality assumption? What about variables/measures that are **not** continuous (e.g. counts) and therefore can't be normal? `Monotonicity` asks do the two variables tend to increase together, or do $Y$ tend to decrease as $X$ increases.

In the parametric/bivariate normal/linearity context:

$$
\begin{aligned}
    +1 &= \text{perfect positive linear} \\\\
    -1 &= \text{perfect negative linear}
\end{aligned}
$$


In the nonparametric/monotonicity settings, we'd like:

$$
\begin{aligned}
    +1 &= \text{perfect increasing monotone} \\\\
    -1 &= \text{perfect decreasing monotone}
\end{aligned}
$$


## Correlation in bivariate data
The key idea (again) is **ranks**, and it requires a notion of ordering. **Exact tests** are based on simulation of the permutation type. A simple scheme would be two paired samples (measurements) with $n$ observations on each measurement. If we fix the order of one of the variables:

| V1 (ranks) | V2   |
| ---------- | ---- |
| 1          | ?    |
| 2          | ?    |
| $\vdots$   | ?    |
| $n$        | ?    |

 and look at all possible orderings of the ranks for the second variable - there are $n!$ of them. Compute measure of correlation for each of these to build the empirical distribution. 

### Spearman rank correlation coefficient
A popular measure is the `Spearman rank correlation coefficient`. It's essentially `Pearson's correlation` calculated on the **ranks** instead of the raw data. Some notations:

- $\rho_s$ - population value
- $r_s$ - sample value
- $(x_i, y_i)$ - paired observations, $i = 1 \cdots n$
- $r_i$ - ranks assigned to $x$ values $i = 1 \cdots n$
- $s_i$ - ranks assigned to $y$ values $i = 1 \cdots n$

$$
r_s = \frac{\sum\limits_{i=1}^n{(r_i - \bar{r})(s_i - \bar{s})}}{\sqrt{\sum\limits_{i=1}^n{(r_i - \bar{r})^2} \sum\limits_{i=1}^n {(s_i - \bar{s})^2}}}
$$



If there are no ties,

$$
r_s = 1 - \frac{6T}{n(n^2-1)} \quad T = \sum\limits_{i=1}^n(r_i - s_i)^2
$$


If $r_i = s_i$ for all $i$, i.e. ranks on $x$ are equal to the ranks on $y$, $T = 0$ and $r_s = 1$. If ranks are perfectly reversed:

| $x$  | 1    | 2     | $\cdots$ | $n$  |
| ---- | ---- | ----- | -------- | ---- |
| $y$  | $n$  | $n-1$ | $\cdots$ | 1    |

We have a perfect monotonically decreasing trend. Here $r_i + s_i = n + 1$ for all $i$. We can show that in this case $r_s = -1$. Intermediate cases (not perfect monotone decreasing/increasing) $\Rightarrow r_s$ is somewhere between $-1$ and $+1$.

### Kendall rank correlation coefficient
The other widely used measure is `Kendall's tau` ($\tau$), which is built on the `Mann-Whitney formulation` (hence also the `J-T test` for ordered alternatives). This is often used as a measure of agreement between judges (how well do two judges agree on their rankings).

We first order the values of the first variable to get $r_i = i$ for all $i$. If there's a positive rank association, ranks for the second variable, $s_i$, should also show an increasing trend; if there's a negative rank association, $s_i$ should show a decreasing trend. 

Order of the $1^{st}$ variable is fixed, $r_i = i$. For the $s_i$, count `concordances` $n_c$  and `discordances` $n_d$, which are pairs that **follow** the ordering and that **reverse** the ordering, respectively. That is, for $i = 1 \cdots n-1$ and $j> i$, count as a `concordance` ($+1$) if $s_j - s_i > 0$ and a `discordance` ($-1$) if $s_j - S-i < 0$. Our test statistic is

$$
t_k = \frac{n_c - n_d}{n(n-1)/2}
$$


**Example**: We have scores on two exam questions for 12 students:

| $Q_1(x)$     | $r_i$ | $Q_2(y)$     | $s_i$ |
| ------------ | ----- | ------------ | ----- |
| 1            | 1     | 13           | 2     |
| 3            | 2     | 15           | 3     |
| 4            | 3     | 18           | 5     |
| 5            | 4     | 16           | 4     |
| 6            | 5     | 23           | 6     |
| 8            | 6     | 31           | 7     |
| 10           | 7     | 39           | 9     |
| 11           | 8     | 56           | 12    |
| 13           | 9     | 45           | 11    |
| 14           | 10    | 43           | 10    |
| 16           | 11    | 37           | 8     |
| 17           | 12    | 0            | 1     |
| Scores 12/20 |       | Scores 12/60 |       |

For $s_1 = 2$, we can get the number of concordant pairs by counting how many $s_j (j > 1)$ are bigger than $2$. There are $10$ concordant pairs and $1$ discordant pair. Similarly, for $s_3 = 5$ there are $7$ concordant pairs and $2$ discordant pairs.

In total,
$$
\begin{gather*}
    n_c = 10 + 9 + 7 + 7 + 6 + 5 + 3 + 0 + 0 + 0 + 0 = 47 \\\\
    n_d = 1 + 1 + 2 + 1 + 1 + 1 + 2 + 4 + 3 + 2 + 1 = 19
\end{gather*}
$$


Note that we have $n_c + n_d = \frac{1}{2}n(n-1)$. If all pairs are concordant, $n_d = 0$ and $t_k = 1$; if all pairs are discordant, $n_c = 0$ and $t_k = -1$. If there is a mix of concordant and discordant pairs, $t_k$ will range between $-1$ and $+1$. In the example, $t_k = \frac{47 - 19}{66} = 0.4242$. 

### Asymptotic results
The number of concordances is essentially the `J-T statistic`. In large enough $n$ (opinion differs on how large is large), we can use an asymptotic normal distribution instead of the exact permutation-based distributions. There haven been numerous approaches suggested in the literature:

$$
\begin{aligned}
    \text{Spearman} && z &= r_s \sqrt{n-1} \dot{\sim} N(0, 1) && \text{under } H_0: \rho_s = 0 \\\\
    && t &= \frac{r_s \sqrt{n-2}}{\sqrt{1 - r_s^2}} \dot{\sim} t_{n-2} && \text{under } H_0: \rho_s = 0 \\\\
    \text{Kendall} && z &= \frac{3t_k\sqrt{n(n-1)}}{\sqrt{2(2n+5)}} \dot{\sim} N(0, 1) && \text{under } H_0: \tau_k = 0
\end{aligned}
$$


**A few remarks:**

1. These approximations are for $H_0$: correlation is 0. Testing other values of potential interest is harder in general. A simulation based approach may help.
2. Confidence intervals for the population correlation are also "tricky" using asymptotic results, because the limits may fall outside the $[-1,1]$ interval. 
3. If the data are really from a bivariate normal distribution `BVN`, both $r_s$ and $t_k$ have pretty high efficiency relative to the Pearson coefficient - $0.912$ in both cases. {{<hl>}}The Kendall coefficient tends to be more powerful than Pearson when the data are from long-tailed, symmetric distributions.{{</hl>}} 

## Ranked data for several variables
We can also consider ranked data for several variables, e.g. more than two judges, and see how well they agree with each other. We often want to test for evidence of **concordance** between rankings of the units. Concordance, whether measured by Kendall's statistic or the Friedman modification, will be one sided: if $H_0$ of no association is rejected, we'll be in the direction of positive association.

### Kendall's W
Kendall's W is a normalization of the statistic of the [Friedman test]({{< ref "/courses/nonparam-stat/3-multiple-samples/3.3-three-or-more-samples/index.md#the-friedman-test" >}}). Let's consider the following example for concordance for multiple judges. Suppose four judges rank five items as follows:

| $_\text{Item}\backslash ^\text{Judge}$ | A    | B    | C    | D    | Sum  |
| -------------------------------------- | ---- | ---- | ---- | ---- | ---- |
| I                                      | 1    | 5    | 1    | 5    | 12   |
| II                                     | 2    | 4    | 2    | 4    | 12   |
| III                                    | 3    | 3    | 3    | 3    | 12   |
| IV                                     | 4    | 2    | 4    | 2    | 12   |
| V                                      | 5    | 1    | 5    | 1    | 12   |

Judges A and C have perfect concordance with each other. Same for judges B and D. However, $(A, C)$ and $(B, D)$ have perfect discordance, and all items end up with the **same** total rank sum. A measure based directly on, say, the Friedman test (or anything that just considers the total rank sums) can't detect a pattern like this.

The test statistic is based on comparison of the rank sums for each unit, as in the Friedman test, but **scaled** by the maximum value attainable (which happens when the judges are in perfect agreement with each other). This maximum value depends both on the number of items and the number of judges:

$$
W = \frac{S}{\max(S)}
$$


where $S$ is the sum of squares of deviations of rank sums from mean under random allocation. The scaling $\max(S)$ allows us to calibrate and interpret the statistic in a meaningful way. $W$ takes value $1$ when there is complete agreement (because $S = \max(S)$ by definition), and takes value $0$ when there is no agreement **or** contrary opinions are held by pairs of judges (as in example above). This is `Kendall's W`, also known as Kendall's coefficient of concordance.

If there are no ties, $\max(S) = \frac{m^2n(n^2 - 1)}{12}$ where $m$ is the number of judges and $n$ is the number of items being judged.

### Another example
There are six contestants in a diving competition. Three judges each independently ranked their performances in order of merit ($1$ for best, $6$ for worst):

| $_\text{Contestant}\backslash ^\text{Judge}$ | A    | B    | C    | Rank total |
| -------------------------------------------- | ---- | ---- | ---- | ---------- |
| I                                            | 2    | 2    | 4    | 8          |
| II                                           | 4    | 3    | 3    | 10         |
| III                                          | 1    | 1    | 2    | 4          |
| IV                                           | 6    | 5    | 5    | 16         |
| V                                            | 3    | 6    | 1    | 10         |
| VI                                           | 5    | 4    | 6    | 15         |

The average rank sum in this case is $10.5$:

$$
S = (8-10.5)^2 + (10-10.5)^2 + \cdots + (15-10.5)^2 = 99.5
$$


To interpret this value, we need to scale it by $\max(S)$:

$$
\begin{gather*}
    \max(S) = \frac{3^2 \times 6(6^2 - 1)}{12} = 157.5 \\\\
    W = \frac{99.5}{157.5} = 0.632
\end{gather*}
$$


Which is a fair amount of agreement. If we want a p-value for $W$, there are some tables for it. We can also use the `Friedman test` for p-value calibration. We can show that

$$
W = \frac{T}{m(n-1)}
$$


where $T$ is the Friedman test statistic. In this example, the p-value is $0.062$ for $H_0$: no concordance, so we may conclude there is agreement though we're unsure how solid it is.

Another possibility for an overall measure of agreement is to look at **all pairwise** measures of correlation/agreement. With $m$ judges, there are $\frac{1}{2}m(m-1)$ pairs. With Spearman's rank correlation coefficient for instance, if we take the mean of the Spearman rank correlations across the pairs,

$$
R_s = \frac{mW - 1}{m-1}
$$


To interpret / calibrate this, we can think of the extreme cases. If all pairwise rankings are in perfect agreement, the pairwise $r_s$ are all 1, which makes $R_s = 1$ and thus $W = 1$. When $W=0$ (no agreement / discordance among pairs), $R_s = \frac{-1}{m-1}$, which can be shown to be the least possible value of $R_s$.