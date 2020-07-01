---
title: "Methods for Paired Samples"
slug: "nonparametric-methods-paired-samples"
categories:
  - Nonparametric Methods
tags:
  - Nonparametric Methods
  - Statistics
  - Paired Samples
summary: "An obvious extension of the one-sample procedures."
date: 2019-04-29T14:22:47-04:00

toc: true  # Show table of contents? true/false

weight: 31
---

We now proceed to an obvious extension of the one-sample procedures: paired samples. This chapter is going to be much shorter as the procedures are very similar to the previously discussed one-sample tests.

Paired observations are pretty common in applications as each unit acts as its own control. For instance, before treatment vs. after, left side vs. right side, identical twin A vs. identical twin B, etc. Another way this comes about is via `matching` (a common experimental design technique): the researcher matches each unit with another unit that shares relevant (or maybe irrelevant) characteristics - age, education, height, etc. and the two units differ on the treatment. Then, use (some of) our one-sample methods on the **differences** between the paired observation values. By analogy, consider the paired t-test, which looks very much like a one-sample t-test.

## Wilcoxon signed rank test
Suppose 12 sets of identical twins were given a psychological test to measure the amount of aggressiveness in each person's personality. We're interested in whether the firstborn tends to be more aggressive than the other. A higher score is indicative of higher aggressiveness. We hypothesize $H_0$: the firstborn does not tend to be more aggressive (no difference is also okay), versus $H_1$: the firstborn twin tends to be more aggressive.

Here's the data we need for the signed rank test on the differences $d_i$:

| Firstborn ($x_i)$ | Secondborn ($y_i$) | Difference ($d_i$) | Rank | Signed rank |
| ----------------- | ------------------ | ------------------ | ---- | ----------- |
| 86                | 88                 | 2                  | 4    | 4           |
| 71                | 77                 | 6                  | 8    | 8           |
| 77                | 76                 | -1                 | 2.5  | -2.5        |
| 68                | 64                 | -4                 | 5    | -5          |
| 91                | 96                 | 5                  | 6.5  | 6.5         |
| 72                | 72                 | 0                  | (1)  | $-$         |
| 77                | 65                 | -12                | 11   | -11         |
| 91                | 90                 | -1                 | 2.5  | -2.5        |
| 70                | 65                 | -5                 | 6.5  | -6.5        |
| 71                | 80                 | 9                  | 10   | 10          |
| 88                | 81                 | -7                 | 9    | -9          |
| 87                | 72                 | -15                | 12   | -12         |

We use mid-ranks for the ties, then use [Sprent and Smeeton "device"]({{< relref "../../2-single-samples/2.1-loc-inf-single-sample/index.md#deviations-of-zero" >}}) for handling the deviation of 0. Our test statistic:

$$
S_+ = 4 + 8 + 6.5 + 10 = 28.5
$$


Now we can get the normal approximation of the test statistic:

$$
\frac{S_+ - \frac{1}{4}n(n+1)}{\sqrt{\frac{1}{24}n(n+1)(2n+1)}} = \frac{28.5 - \frac{1}{4} \times 12 \times 13}{\sqrt{\frac{1}{24} \times 12 \times 13 \times 25}} = -0.824
$$


And the score representation:
$$
\frac{S_+ - \frac{1}{2} \sum\limits_{i=1}^n {|S_i|}}{\frac{1}{2}\sqrt{\sum\limits_{i=1}^n {S_i^2}}} = \frac{28.5 - \frac{1}{2} \times 77}{\frac{1}{2}\sqrt{648}} = -0.786
$$


For the score representation, we see a slight difference when we compare to: $\sum\limits_{i=1}^{12} i = 78$, and $\sum\limits_{i=1}^{12} i^2 = 650$. The difference is caused by the $-$ rank in the case where $d_i = 0$. In either case, there's really no evidence that the firstborn twin is more aggressive. 

### Remarks
1. Within a pair, observations may not be independent (obviously, as almost by construction sometimes they won't be), but the **pairs** themselves should be.
2. The typical $H_0$ here would be that the **median of the differences** is $0$. If the differences are assumed to have a **symmetric distribution** (about $0$), then the Wilcoxon approach is suitable. Note then that the individual unit measurements do **not** need to be assumed symmetric. You need to think about whether this is realistic for your given situation.
3. The alternative therefore refers to a **shift in centrality** of the differences. You need to think about whether this is an interesting, relevant, or important question.

## McNemar's test
The test is a less obvious use of, or a modification to the sign test. Here's example 5.4 in the book. We have records on all attempts of two rock climbs, successful or not. For the $108$ who tried both:

|                      | First climb success | First climb failure |
| -------------------- | ------------------- | ------------------- |
| Second climb success | 73                  | 14*                 |
| Second climb failure | 9*                  | 12                  |

Is there evidence that one climb is harder? The only ones that had information **for this question** are the ones that succeeded in one but failed in another. The people who succeeded / failed in both are essentially "ties".

We can frame this as a Binomial setup. Think of success $\rightarrow$ failure as a "+" for the first climb, and failure $\rightarrow$ success as a "-" for the first climb. This puts us under a sign test situation. Under $H_0:$ climbs are equally difficult, and the probability of a "+" is the same as the probability of a "-" $\Rightarrow Bin(23, 0.5)$ with a "+". The p-value is 0.405, so there seems to be no difference in the difficulty of the climbs.

This version of the sign test is called `McNemar's test`. In general, we have pairs of data $(x_i, y_i)$, $i=1 \cdots m$ where $X_i$ and $Y_i$ take values 0 and 1 only. There are 4 patterns of outcomes: $(0, 0), (0, 1), (1, 0), (1, 1)$. The paired data are summarized in a $2 \times 2$ `contingency table` showing the classifications of $Y_i$ and $X_i$:

|           | $Y_i = 0$                       | $Y_i = 1$                       |
| --------- | ------------------------------- | ------------------------------- |
| $X_i = 0$ | $a\text{ (\# of (0, 0) pairs)}$ | $b\text{ (\# of (0, 1) pairs)}$ |
| $X_i = 1$ | $c\text{ (\# of (1, 0) pairs)}$ | $d\text{ (\# of (1, 1) pairs)}$ |

### Assumptions
1. Pairs are mutually independent. 
2. Two categories for each outcome ($X_i = 0, 1\,/\,Y_i = 0, 1$)
3. The difference $P(X_i = 0, Y_i = 1) - P(X_i = 1, Y_i = 0)$ is negative for all $i$, 0 for all $i$, or positive for all $i$.

The null hypothesis can be formulated in various equivalent ways:

$$
\begin{aligned}
    &H_0: P(X_i = 0, Y_i = 1) = P(X_i = 1, Y_i = 0) && \text{for all } i \\\\
    \text{vs. } &H_1: P(X_i = 0, Y_i = 1) \neq P(X_i = 1, Y_i = 0) && \text{for all } i
\end{aligned}
$$


or

$$
\begin{aligned}
    &H_0: P(X_i = 0) = P(Y_i = 0) && \text{for all } i \\\\  \text{vs. } &H_1: P(X_i = 0) \neq P(Y_i = 0) && \text{for all } i
\end{aligned}
$$


or 

$$
\begin{aligned}
    &H_0: P(X_i = 1) = P(Y_i = 1) && \text{for all } i \\\\
    \text{vs. } &H_1: P(X_i = 1) \neq P(Y_i = 1) && \text{for all } i
\end{aligned}
$$


The first one, in terms of our example, is saying $P(\text{success} \rightarrow \text{failure})$ is equal to $P(\text{failure} \rightarrow \text{success})$. The second one is saying $P(\text{failure on first climb})$ $=$ $P(\text{failure on second climb})$.

Our test statistic - if $b$ and $c$ are reasonably small, take $T_2 = b$. Under $H_0$, $b \sim Bin(b+c, 0.5)$, which is our usual sign test. Or for larger $b, c$:

$$
T_1 = \frac{(b-c)^2}{b+c}
$$


Note that it does **not** depend on $a$ and $d$. Under $H_0$, $T_1$ is approximately $\chi^2_{(1)}$. To derive $T_1$, let's call $b + c = n$. If $n$ is large enough, we can use a normal approximation to the Binomial:
$$
\frac{T_2 - \frac{1}{2}n}{\sqrt{\frac{1}{2} \cdot \frac{1}{2}n}} = \frac{b - \frac{1}{2}(b+c)}{\frac{1}{4}(b+c)} = \frac{b-c}{\sqrt{b+c}} = \mathbf{Z} \approx N(0, 1)
$$

$$
T_1 = \mathbf{Z}^2 \Rightarrow T_1 \approx \chi^2_{(1)}
$$


And that's it! Paired-sample tests really aren't that different from one-sample tests. Now we're ready for some real problems.
