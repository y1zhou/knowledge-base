---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Fundamentals of Nonparametric Methods"
slug: "nonparametric-methods-fundamentals"
summary: "Some basic tools such as the permutation test and the binomial test. We also introduce order statistics and ranks, which will come in handy in later chapters."
categories:
  - Nonparametric Methods
tags:
  - Nonparametric Methods
  - Statistics
date: 2019-01-25T22:50:34-05:00
lastmod: 2019-01-25T22:50:34-05:00
draft: false  # Is this a draft? true/false
toc: true # Show table of contents? true/false

featured: true
weight: 12

header_image:
  unsplash_id: "NL_DF0Klepc"
---

Here we continue to present some basic tools (permutation test and sign test) and principles (order statistics, ranks, and efficiency) that will be useful moving forward.

## Permutation (Randomization) test

Consider the following pretty artificial scenario in which four out of nine subjects are selected **at random** to receive a new drug, and the other five get a placebo. After some time, all nine subjects are assessed on some outcome and are **ranked** from "best" ($rank = 1$) to "worst" ($rank = 9$). We assume no ties - for now. 

Here's our question: if the new drug has no beneficial effect, what is the probability that the subjects who got it were ranked $(1, 2, 3, 4)$, i.e. best responses?

First, we need to think about what it means when we say that the four were selected "at random". If the drug really has no effect, it means that the labels "new drug" and "placebo" are essentially meaningless, and **any** subject is equally likely to be ranked low, medium, or high after the treatment. The ranks are essentially assigned at random.

There are in total $\binom{9}{4} = 126$ ways of picking four people out of the nine to receive the new drug. If the drug has no effect, then the set of ranks belonging to the chosen four is **equally likely** to be any of the 126 possible sets of ranks. Here are some of the possibilities:

| A        | B        | C        | D        |                                               |
| -------- | -------- | -------- | -------- | --------------------------------------------- |
| 1        | 2        | 3        | 4        | $\rightarrow$ most favorable to the new drug  |
| 1        | 2        | 3        | 5        |                                               |
| 1        | 2        | 3        | 6        |                                               |
| 1        | 2        | 3        | 7        |                                               |
| $\vdots$ | $\vdots$ | $\vdots$ | $\vdots$ |                                               |
| 5        | 6        | 7        | 9        |                                               |
| 5        | 6        | 8        | 9        |                                               |
| 5        | 7        | 8        | 9        |                                               |
| 6        | 7        | 8        | 9        | $\rightarrow$ least favorable to the new drug |

We could think of testing

$$
\begin{aligned}
    & H_0: \text{new drug has no effect} \\\\
    vs. \\, &H_1: \text{new drug has an effect}
\end{aligned}
$$

using the enumeration of the rank outcomes. This is the basic logic of the `permutation test`. If there is no difference / no effect ($H_0$), then the labels are arbitrary and the ranks likewise can be thought of as arbitrarily assigned, so we can look at the configuration of ranks actually observed and see how likely it is. With the two-sided alternative $H_1$, we consider both extremes $(1, 2, 3, 4)$ and $(6, 7, 8, 9)$, yielding a probability of $\frac{2}{126} \approx 0.016$. 

Usually, instead of working with the ranks, we work with a **function of the ranks** that has a low value if all of the ranks are low, a high value if all of the ranks are high, and an intermediate value if there is a mix of high and low ranks (in the "treated group"). For instance, **sum (or average) of the ranks**. For our example, the sum will range from 10 to 30:

| $S \text{ (rank sum)}$ | # of occurrences | Probability |
| ---------------------- | ---------------- | ----------- |
| 10                     | 1                | 0.008       |
| 11                     | 1                | 0.008       |
| 12                     | 2                | 0.0159      |
| 13                     | 3                |             |
| 14                     | 5                |             |
| $\vdots$               | $\vdots$         | $\vdots$    |
| 26                     | 5                |             |
| 27                     | 3                |             |
| 28                     | 2                |             |
| 29                     | 1                |             |
| 30                     | 1                |             |

The distribution is symmetric. We could also define a `rejection region`, as in classical testing, by looking for values of $S$ in each tail with $prob. < 0.025$ (for two-tailed test, and taking $\alpha = 0.05$), to get an overall level of $0.05$. This gives $S = \{10, 11, 29, 30\}$ with cumulative prob. $\approx 0.0317 > 0.025$.

> 1. We can't attain the exact 0.05 level here. This is characteristic of permutation tests because we're dealing with discrete distributions for the test statistic. As the sample sizes increase, you can get closer.
> 2. We made **no assumptions** about what distribution the outcome was drawn from. But note that we *also* don't use the actual outcome values themselves - just the ranks! There is a loss of information for a gain in flexibility and robustness (**tradeoff in this case**).

## Binomial tests
Another key component is tests that are built on the `binomial distribution`. We'll introduce this with an example.

### Motivating example
We have data on survival times of ten patients with a certain type of cancer. But, for one patient the precise survival time is not known - the study only followed subjects for 362 weeks and that subject is still alive at that point.

In this case, we say that the observation is **censored**. The data are (survival time in weeks):

\\[49, 58, 75, 110, 112, 132, 151, 276, 281, 362^*\\]

> $*$: censored, survival time > 362 weeks

Suppose we want to test the hypothesis that the **median survival time** (in the relevant population) is 200 weeks, vs. the alternative that it is not. Let $\theta$ be the median survival time in the population from which the sample was drawn. We're interested in:

$$
\begin{aligned}
    &H_0: \theta = 200 \\\\
    vs. \\, &H_1: \theta \neq 200
\end{aligned}
$$

For classical (parametric) approaches, this scenario has two main complications:

1. the censored observation need special methods, and
2. the fact that we're looking at the median rather than the mean.

The `Binomial test` offers a way around *both* of these issues simultaneously. If we have a random sample from any continuous distribution with median 200, each sample value is equally likely to be above or below 200. Define a "success" to be a value above 200, and a "failure" to be a value below 200. Under $H_0$, a success is just as likely as a failure, since the median $\theta$ is 200 in the **population**.

So, in a sample of size 10, the number of successes is $Bin(10, 0.5)$. We have 3 successes, including the censored observation (whose precise value no longer matters!).
$$
\text{F, F, F, F, F, F, F, S, S, S}
$$
Under $H_0$, $P \text{(3 successes)} = {10 \choose 3} 0.5^{10}$. We also need to consider the "more extreme" outcomes in this tail: 2 or 1 or 0 successes. In total:
$$
{10 \choose 3}0.5 ^{10} + {10 \choose 2}0.5 ^{10} + {10 \choose 1}0.5 ^{10} + {10 \choose 0}0.5 ^{10} \approx 0.1719
$$


We have a two-sided test, so we also need to consider departures from $H_0$ in the other tail. Since $Pr(\text{success}) = Pr(\text{failure}) = 0.5$, this distribution is symmetric, so the p-value is $2 \times 0.1719 = 0.3438$.

This Binomial test is also called the `sign test`, and is very popular as a nonparametric test due to its simplicity.

### Intepreting results
As in the previous example of the permutation test, because of the discreteness of the binomial, not all levels (p-values) are attainable. As $n$ increases, this becomes less of an issue, but for small or moderate sample sizes, this argues against **strict** cutoffs (fixed, but essentially arbitrary) such as "p-value $< 0.05$ $\Rightarrow$ 'statistically significant'".

It is better to report the p-value itself, and if possible, also report a `confidence interval` about the population quantity of interest. With $(0, 1, 9, 10)$ successes (low and high # of successes), we have a total probability of 0.0216​, which is "reasonable evidence" against $H_0$. If we add in 2 and 8, that probability jumps to ​0.1094, which you may or may not be willing to take as evidence against $H_0$.

If we go with 0.0216, we have between 2 and 8 successes if the median specified in $H_0$ has any value greater than 58 but less than 281. So $(58,281)$ is a $97.84\\%$ confidence interval for $\theta$, the population median. This is a **very** wide `confidence interval` - probably too wide for a clinical setting. We've gained flexibility in analysis, but we've thrown out a lot of the information in the data.

### Asymmetric cases
In a dental practice, experience has shown that $75\\%$ of adult patients require treatment after a routine checkup. So, in a sample if 10 independent patients, the number, $S$, needing additional treatment is $Bin(10, 0.75)$. The probabilities for each outcome are:

| $S\text{(\# of successes)}$ | prob. of individual outcomes |
| --------------------------- | ---------------------------- |
| 0                           | 0.0000                       |
| 1                           | 0.0000                       |
| 2                           | 0.0004                       |
| 3                           | 0.0031                       |
| 4                           | 0.0162                       |
| 5                           | 0.0584                       |
| 6                           | 0.146                        |
| 7                           | 0.2503                       |
| 8                           | 0.2816                       |
| 9                           | 0.1877                       |
| 10                          | 0.0563                       |

We can see that the distribution is no longer symmetric in the two tails! Suppose we had a sample of size 10 from a *different* dental practice and wanted to test: 


$$
H_0: p= 0.75 \text{ vs. } H_1: p > 0.75
$$


The most extreme[^1] result would be 10 successes in the second practice, which under $H_0$ ($p = 0.75$) has probability 0.0563, which is "just above" the "standard" $\alpha = 0.05$ threshold. So with a traditional testing approach, we could never reject $H_0$. Another situation where it makes more sense to report the p-value itself. The other one-tailed test:


$$
H_0: p= 0.75 \text{ vs. } H_1: p < 0.75
$$


Here we keep accumulating the lower tail[^2] probabilities until we hit the p-value threshold. In this case, we can reject $H_0$ in the traditional framework.

When the distribution is **asymmetric** and we have a two-tailed test, there are two options:
1. Find the point in the other tail with equal or lower cumulative probability, e.g. both tail probabilities should approach 0.025 to get a p-value of 0.05.
2. Take tails eqidistant from the mean (in discrete cases, take the same number of bins).

## Order statistics and ranks

Ranking is at the basis of many nonparametric tests. It's a major way of dropping distributional assumptions at the cost of losing some of the information in the raw data. We can use ranks / ordered data.

In general, if we have observations $x_1, x_2, \cdots, x_n$ from a continuous distribution (i.e. no ties), we denote the `order statistics` $x_{(1)} < x_{(2)} < \cdots < x_{(n)}$:

- $x_{(1)}$ is the smallest obs. in the sample.
- $x_{(1)} < x_{(2)} < \cdots < x_{(n)}$ is the sample ordered from smallest to largest.

The **median** of the sample can be defined in terms of the order statistics:


$$
\begin{aligned}
    &\text{n is odd} \Rightarrow \text{sample median is }x_{\left(\frac{n+1}{2}\right)} \\\\
    &n = 2m \Rightarrow \text{sample median is } \frac{1}{2}\left(x_{(m)} + x_{(m+1)}\right)
\end{aligned}
$$
We can also define `measures of dispersion` - range or interquartile range - in terms of the order statistics, e.g. $x_{(n)} - x_{(1)}$ is a possible measure of dispersion. Another key use of the order statistics is to build the `empirical distribution function`:


$$
S_n(x) = \frac{1}{n} (\text{\# sample values} \leq x)
$$
In terms of the order statistics, 
$$
S_n(x) = \begin{cases}
    0, & x < x_{(1)}  \\\\
    \frac{i}{n}, & x_{(i)} \leq x < x_{(i+1)}, & i = 1 \cdots n-1 \\\\
    1, & x \geq x_{(n)}
\end{cases}
$$
which is a "step function" with jumps of $\frac{1}{n}$ at each of the observed data points. This is an estimator of the population CDF $F(x)$.

## Power and efficiency
How effective is a procedure in using the information in the sample? {{< hl >}}In general, nonparametric methods are less efficient than parametric counterparts when the assumptions of the parametric approaches are met.{{</hl >}} Nonparametric methods, e.g. based on ranks, replace the actual observed data with ordered values $\Rightarrow$ toss out information $\Rightarrow$ lose power / efficiency.

One way of looking at efficiency is through `asymptotic relative efficiency (ARE)`. "Asymptotic" meaning as sample sizes increase, and "relative efficiency" is comparing two procedures.

Consider two sequences of tests, $T_1$ and $T_2$ (different tests with increasing sample sizes, e.g. 5, 6, 7, ... samples) where $\alpha$ (probability of [Type I error](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors#Type_I_error)) is fixed. We let $H_1$ vary in a way that $\beta$ (probability of [Type II error](https://en.wikipedia.org/wiki/Type_I_and_type_II_errors#Type_II_error)) remains constant as the sample size for test sequence $T_1$, call it $n_1$, increases. For each value of $n_1$, the idea is to determine $n_2$ such that the test sequence $T_2$ has the same $\beta$ for the particular alternative.

The proposed idea is to test what are the sample sizes you'd need for tests $T_1\,(n_1)$ and $T_2\,(n_2)$ to get the **same test performance**. The closer $n_1$ and $n_2$ are to each other, the tests have closer efficiency - they use the data with (near) equal effectiveness.

A bit more rigorously, a bigger sample size usually leads to increased power for alternatives closer to the null. With a larger sample, you can detect smaller differences. For big samples, the ratio $\frac{n_1}{n_2}$ is potentially informative, and it can be shown that (under some circumstances) $\frac{n_1}{n_2}$ tends to a limit as $n_1 \rightarrow \infty$. This is the `asymptotic relative efficiency`. Crucially, nonparametric methods *can be* more powerful than their parametric counterparts when the assumptions of the latter don't hold.

Next, we'll discuss [location inference on a single sample]({{< ref "/series/nonparam-stat/2-single-samples/2.1-loc-inf-single-sample/index.md" >}}), and the tool we'll be using is the Wilcoxon signed-rank test.

[^1]: One-sided alternative values in the upper tail (high # of success) $\Rightarrow$ favorable to $H_1$.
[^2]: One-sided alternative values in the lower tail (low # of successes) $\Rightarrow$ favorable to $H_1$.
