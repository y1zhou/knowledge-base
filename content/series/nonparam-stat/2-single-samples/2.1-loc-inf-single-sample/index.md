---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Location Inference for Single Samples"
slug: "nonparametric-methods-single-sample-location-inference"
categories:
  - Nonparametric Methods
tags:
  - Nonparametric Methods
  - Statistics
summary: The Wilcoxin signed rank test explained.
date: 2019-03-26T21:12:45-05:00
lastmod: 2019-03-26T21:12:45-05:00
draft: false  # Is this a draft? true/false
toc: true  # Show table of contents? true/false

weight: 21
---

Previously we've used the sign test to look at the median (a measure of location) survival time with a censored data point. The original observations were transformed into "successes" and "failures" and a **lot** of information was thrown out.

For location inference, a basic one-sample procedure is the `Wilcoxon Signed-Rank test`. It's used to test whether a sample comes from a population with a specified mean or median.

## Wilcoxon signed-rank test
Suppose $n$ observations $x_1, x_2, \cdots, x_n$, are a sample from a **symmetric continuous** distribution with unknown median, $\theta$. We want to test:

$$
H_0: \theta = \theta_0 \text{ vs. } H_1: \theta \neq \theta_0
$$

The assumption of symmetry implies that under $H_0$:

$$
\begin{aligned} d_i = x_i - \theta_0 && i = 1, 2, \cdots, n \end{aligned}
$$

- Each $d_i = x_i - \theta_0$ is equally likely to be positive or negative (i.e. each $x_i$ is equally likely to be above or below $\theta_0$ under $H_0$ - logic for sign test).
- The magnitude $|d_i|$ of any size is equally likely to be positive or negative.

A symmetric population distribution has a mean that coincides with the median. In those circumstances the test may also be formulated in terms of means. The test procedure is simple: 
1. Calculate the discrepancies of each observation from $\theta_0$, the median hypothesized under $H_0$. 
2. Order these magnitudes (i.e. in absolute value) from smallest ($\text{rank} = 1$) to largest ($\text{rank} = n$). 
3. Assign a $+$ sign to ranks corresponding to $d_i > 0$, and a $-$ sign to ranks corresponding to $d_i < 0$.

From here, we can define a few different test statistics. We denote by $S_+$the the sum of positive ranks, and by $S_-$ the sum of ranks associated with negative deviations. These two are equivalent to each other because the total sum of ranks in sample of size $n$ is fixed: $1 + 2 + \cdots + n = S_+ + S_-$. We may use any of the statistics $S_+, S_-$, or a third statistics $S_d = |S_+ - S_-|$ as a test statistic.  All three have the same information about the plausibility of $H_0$. Under $H_0$, we'd expect $S_+$ and $S_-$ to be roughly equal, and likewise $S_d$ should be close to $0$. Intermediate values of $S_+$ or $S_-$ are more supportive of $H_0$.

We can use the permutation test approach to get the **exact** distribution for any of the three test statistics. We look at all the possible allocations of $+/-$ signs to the ranks $1,2, \cdots, n$. Let's take a look at the following example.

## Heart rate example
Heart rate (beats per minute) when standing was recorded for seven people. Assume a symmetric distribution for heart rate in the population. Continuity is questionable as heart rate is an integer, but we're okay if there are no ties, which is usually why we assume continuity.

The observed data is:

$$
73, 82, 87, 68, 106, 60, 97
$$

Suppose we want to test

$$
\begin{aligned}
    &H_0: \theta = 70 \\\\
    \text{vs. } &H_1: \theta > 70 && \text{one-sided alternative}   \end{aligned}
$$

First, compute $|d_i| = |x_i - \theta_0| = |x_i - 70|$, $i = 1, \cdots, 7$:



| Magnitude $\|d_i\|$ | 3    | 12   | 17   | 2    | 36   | 10   | 27   |
| ------------------- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| Rank                | 2    | 4    | 5    | -1   | 7    | -3   | 6    |

Note now we added a minus sign to the ranks of the observations that were smaller than 70. Now we can calculate the test statistics:

$$
\begin{aligned}
    \Rightarrow S_+ &= 2 + 4 + 5 + 7 + 6 = 24 \\\\
    S_- &= 1 + 3 = 4 \\\\
    S_d &= |S_+ - S_-| = 20
\end{aligned}
$$


> Sanity check: $1 + \cdots + 7 = 28 = S_+ + S_-$.

Under $H_0$ ($\theta = 70$), we'd expect $S_+ \approx S_-$ and $S_d \approx 0$. Under $H_1$ ($\theta > 70$), we'd expect $S_+ > S_-$ and $S_d$ to be larger. Here we observe $S_+$ to be appreciably larger than $S_-$, which is in support of $H_1$.

For the permutation test, we need to build the permutation distribution. Note that $S_+$ can take values from 0 (when all ranks are negative $\Rightarrow$ all observations are less than $\theta_0$) up to 28 (all ranks are positive). There are $2^7 = 128$ ways of allocating $+/-$ signs to the ranks. All of them are equally likely under $H_0$. Below we have a few possible configs:

| 1        | 2        | 3        | 4        | 5        | 6        | 7        | $S_+$    | $S_-$    | $S_d$    |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| -        | -        | -        | -        | -        | -        | -        | 0        | 28       | 28       |
| +        | -        | -        | -        | -        | -        | -        | 1        | 27       | 26       |
| $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ |
| +        | +        | +        | +        | +        | +        | +        | 28       | 0        | 28       |

With an one-sided test, the upper tail $H_1 \Rightarrow$ large values of $S_+$ (equivalently small values of $S_-$) are evidence against $H_0$. We want
$$
\begin{aligned}
    P(S_+ \geq 24) &= P(S_+ = 24) + P(S_+ = 25) + \cdots + P(S_+ = 128) \\\\
    &= \frac{7}{128} \approx 0.0547
\end{aligned}
$$


which is pretty unlikely if $H_0$ holds. The evidence against $H_0$ is sufficiently strong to warrant a similar study with a larger group of students.

## Function in R
In R, the function `wilcox.test()` does the exact permutation test if $n < 50$ and there are no ties in data. Otherwise, a normal approximation is used. It should be noted that two options in `wilcox.test` conflict with each other:
- `exact`: compute exact permutation distribution. Ignored when there are ties.
- `correct`: use [continuity correction](#asymptotic-results) based on the normal distribution.
- `exact` overrides `correct`, so when `exact = T` the value of `correct` doesn't matter.

In our case, `wilcox.test(data, "greater", mu = 70)` gives $S_+ = 24$ (R calls it `V`) and a p-value of $0.05469$. We can also do this by hand with the following procedure.

There are $2^{12} = 4096$ configurations in the permutation distribution. The statistic $S_+$ (or $S_-$) can take on 79 possible values - the values $0$ to $1 + 2 + \cdots + 12 = 78$. We may use `dsignrank(z, n)` to build up the permutation distribution for $S_+$ (equivalently $S_-$)[^1], where `z` are the possible values, and `n` is the sample size.

Another useful application is to get a confidence interval (CI) for the population median (mean): `wilcox.test(x,conf.int = T, conf.level = 0.95)`. The test is based on `Walsh averages`.

## Walsh average
Walsh averages are just pairwise averages of the observations. Set $H_0: \theta = \theta_0$ and get the order statistics $x\_{(1)} <x\_{(2)} < \cdots <x\_{(n)}$. Now, if $x_i < \theta_0$, $X_i - \theta_0$ will have a **negative** rank. In fact, the largest negative deviation will belong to $x_{(1)}$, the next largest to $x_{(2)}$, and so on. Likewise, if $x_j > \theta_0$, $x_j - \theta_0$ will have a positive rank, and the largest positive deviation will belong to $x_{(n)}$.

Now, we consider the paired averages of $x_{(1)}$ with each of $x_{(1)}, x_{(2)}, \cdots, x_{(n)}$.

1. Denote by $-p$ the signed rank of the deviation associated with $x_{(1)}$. Then each $\frac{1}{2}(x_{(1)} + x_{(q)})$ is less than $\theta_0$ when the deviation associated with $x_{(q)}$ also has negative rank, or positive rank less than $p$.
2. If $x_{(q)}$ is the smallest observation with a positive rank greater than $p$, then $|x_{(q)} - \theta_0| > |x_{(1)} - \theta_0|$, and then the average of $x_{(1)}$ and $x_{(q)} > \theta_0$. And also any $x_{(r)} > x_{(q)}$.
3. The number of averages involving $x_{(1)}$ that are less than $\theta_0$ equals the negative rank associated with $x_{(1)}$. Likewise, the number of averages involving $x_{(i)}$ less than $\theta_0$ with each of $x_{(i)}, x_{(i+1)}, \cdots, x_{(n)}$ that are less than $\theta_0$ equals the negative rank associated with $x_{(i)}$.
4. The number of `Walsh averages` less than $\theta_0 = S_-$, and the number of Walsh averages greater than $\theta_0 = S_+$.

We can use this correspondence to get point estimate and CI for the population median (mean). To do this:
1. Order the data from smallest to largest.

2. Compute the pairwise Walsh averages and arrange them in a table. There are 78 unique averages. The Walsh averages increase along each row and column.

3. Point estimator of the population median is the **median of the Walsh averages**, which should be the average of the $39^{th}$ and the $40^{th}$ ordered Walsh averages. In our example, they are both $9 \Rightarrow \hat{\theta} = 9$. The median of the Walsh averages is also called the `Hodges-Lehmann estimator`, having been proposed by Hodges and Lehmann[^2]. 

4. To find a $95\\%$ confidence interval using the Walsh averages, we select as end points values of $\theta$ that will just be acceptable if $P = 0.05$. For one tail, $0.025 \times 4096 = 102.4$, so the cumulative sum can't be more than 102.4. CI puts roughly 2.5% of density in each tail, so we look for $14^{th}$ smallest and $14^{th}$ largest Walsh averages $\Rightarrow (2.5, 16.0)$.

## Asymptotic results
What happens as sample size $n$ increases? The number of possible configurations (possible assignments of $+/-$ signs to the ranks) becomes large very quickly - $2^n$ increases rapidly. Therefore, it's not feasible to get the exact distribution even with the aid of software.

But we also have symmetry, and as $n$ increases, the number of possible values of $S_+$ also increases ($\frac{1}{2}n(n+1)$). There is a **normal approximation** to the exact distribution of the test statistic! We can easily show that 
$$
\begin{gather*}
    E(S_+) = \frac{1}{4}n(n+1) \\\\
    Var(S_+) = \frac{1}{24}n(n+1)(2n+1)
\end{gather*}
$$


since the sum of the integers from $1$ to $n$ is $\frac{n(n+1)}{2}$ and the sum of squares is $\frac{n(n+1)(2n+1)}{6}$. For large enough $n$, we can define a test statistic $Z$


$$
Z = \frac{S_+ - E(S_+)}{\sqrt{Var(S_+)}} \dot\sim N(0,1)
$$


We can also improve the effectiveness of the approximation by using a `continuity correction`, which comes to account for the fact that we approximate a discrete quantity $S_+$ with a continuous distribution (in this case the normal). The idea is if $S$ is the smaller of the sums of the positive of negative ranks, we replace $S$ by $S + \frac{1}{2}$. If $S$ is the larger sum it is replaced by $S - \frac{1}{2}$. In R, this is called with the option `correct = T` in `wilcox.test`.

Another thing we could do when $2^n$ is too large to enumerate all of the possible configurations (e.g. in the millions) is to take a random sample of them instead (e.g. $\approx 10,000$).

## Wilcoxon with ties
We've been assuming the underlying population distribution is continuous, so in theory in the data there should be no ties. In practice, however, ties can of course happen. Real observations never have a distribution that is strictly continuous either because of their nature, or as a result of rounding errors or limited measurement precision. 

The other in theory impossible but in practice often observed are deviations of $0$. We're talking about values in the sample that are exactly equal to the median hypothesis under $H_0$, $\theta_0$. There's no one agreed upon best approach for handling these cases!

### Ties in the observations
One suggestion is to replace the ranks for the tied values with their `mid-ranks`. For example, if we had a sample 


$$
12, 18, 24, 26, 37, 40, 42, 47, 49, 49, 78, 108
$$




We take $H_0: \theta = 30$. The deviations and ranks are then

| $d_i$ | -18  | $-12^*$ | -6   | -4   | 7    | 10   | $12^*$ | 17   | 19$^\dagger$ | 19$^\dagger$ | 48   | 78   |
| ----- | ---- | ------- | ---- | ---- | ---- | ---- | ------ | ---- | ------------ | ------------ | ---- | ---- |
| Rank  | -8   | -5.5    | -2   | -1   | 3    | 4    | 5.5    | 7    | 9.5          | 9.5          | 11   | 12   |



The *s have the same magnitude so their ranks would be tied. They would be ranked $5$ and $6$, so we give them both $5.5$. Similarly the 19's are both ranked 9.5.

We can now calculate $S_+$ and $S_-$ same as before. But the exact distribution changes because of using these mid-ranks. In fact, {{<hl>}}the exact permutation distribution depends on the number of ties and where they fall in the rank sequence.{{</hl>}} It's much harder to work out. The distribution with no ties is unimodal, with values of the statistic confined to ingegral values increasing in steps, while the distribution with ties is heavily multimodal with $S$ taking unevenly spaced and not necessarily integer values. Discontinuities are also more marked.

Streitberg and Rohmel[^3] came up with an algorithm called the `shift algorithm` to handle this situation. The R package `exactRankTests` implements this algorithm under the function `wilcox.exact`. Works for ties or no ties!

Another suggestion is to modify the normal approximation. The idea is to consider each of the signed ranks as a score $S_i$ for observation $x_i$. Under $H_0$, as before each score has equal probability of being positive or negative, so the expected value and the variance of $S_+$ or $S_-$ is 


$$
\begin{aligned}
    E(S_i) &= \frac{1}{2}\sum\limits_{i=1}^{n} |S_i| \\\\
    Var(S_i) &= \frac{1}{4}\sum\limits_{i=1}^n S_i^2
\end{aligned}
$$




These work out to be the same as before for the particular choice of score $\Leftrightarrow$ rank. The **score representation** is thus


$$
Z = \frac{S_+ - \frac{1}{2}\sum{|S_i|}}{\sqrt{\frac{1}{4}\sum{S_i^2}}} \dot\sim N(0,1)
$$


### Deviations of zero
Again, opinions differ for data points that are equal to the median hypothesized under $H_0$. A standard advice is to drop such points from the calculation of $S_+$ (this is equivalent to assigning them rank 0), but this decreases the sample size and we lose data!

Sprent and Smeeton (in part 3.3.4 of the book) proposed a slight alternative: temporarily assign such points a rank of 1 (or the appropriate tied rank if there's more than one zero), then sign-rank all the other observations as usual. Finally, switch the rank(s) associated with zero deviation to 0. This keeps all the data and uses the ranks up to $n$. However, the effect on the exact distribution is unclear.

## Summary
1. The `sign test`, unlike the `Wilcoxon signed rank test`, does not require symmetry of the underlying distribution. When the data come from a skewed distribution (e.g. income), both the `t-test` and the `Wilcoxon` may be inappropriate in the sense that they may not give us valid inference. This depends in part on how skewed the population distribution is. The `sign test` will still be valid. Confidence intervals based on the `t-test` / `Wilcoxon` may also be misleading. Those based on the `sign test` (and `Binomial distribution`) will still be fine. In other situations, all three will lead to similar conclusions.
2. A suggestion: try different analyses and see if your conclusions are consistent.
3. When the symmetry assumption is violated, the `sign test` may have higher efficiency - higher power in tests and shorter confidence intervals for a given confidence level. We can compare the [asymptotic relative efficiency](https://en.wikipedia.org/wiki/Efficiency_(statistics)) of the sign test, Wilcoxon and t-test:
   - ARE of the Wilcoxon compared to the t-test is at least 0.864, and can go up to infinity under some circumstances.
   - The Wilcoxon is never "too bad" and can be very good.
   - When the data are actually normally distributed (the situation where the t-test is optimal), the ARE of Wilcoxon is 0.955, so very little loss here.

Up to this point, we've been talking about location inference - questions about the mean or median of the population distribution, but we can do a lot more. One particularly useful application is studying whether our data are consistent with having been drawn from some specified distribution.



[^1]: For each distribution, there's a series of functions in R to get the density `dsignrank`, distribution function `psignrank`, quantile function `qsignrank` and random numbers `rsignrank`.
[^2]: Hodges Jr, J. L., & Lehmann, E. L. (1963). Estimates of location based on rank tests. *The Annals of Mathematical Statistics*, 598-611.
[^3]: Streitberg, B., & Röhmel, J. (1984). Exact nonparametrics in APL. *ACM SIGAPL APL Quote Quad*, *14*(4), 313-325.
