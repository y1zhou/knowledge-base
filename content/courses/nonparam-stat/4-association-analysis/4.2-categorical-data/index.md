---
title: "Categorical Data"
slug: "nonparametric-methods-categorical-data"
categories:
  - Nonmarametric Methods
tags:
  - Nonmarametric Methods
  - Statistics

summary: "Dealing with contingency tables. Fisher's exact test comes back, together with Chi-squared test and likelihood-ratio test. We also talk about testing goodness-of-fit."
date: 2019-05-06T10:46:18-04:00
toc: true
type: docs  # Do not modify.
weight: 90

menu:
  nonparam-stat:
    name: Categorical Data
    parent: Association Analysis
    weight: 20
---

In this chapter our focus is mostly on "count data" - the data are number of units with particular attributes. The problems appear different at first sight, yet many are solved using procedures developed previously.

`Attributes` are characteristics that the units can have, and are typically qualitative or categorical. Data are arranged in a contingency table. The table has a dimension (row, column, layer, ...) for each attribute. We consider multiple attributes simultaneously. Attributes can be

- Nominal or ordinal
  - `Nominal` - categories are just names, not ordering
  - `Ordinal` - categories can be ordered 
- Explanatory or response variables

The general setup is an $r \times c$ table (rows and columns) where

- $n_{ij}$ is the count in $(i, j)$ cell
- $n_{i+}$ is also denoted $n_{i\cdot}$, and $n_{+j}$ is also denoted $n_{\cdot j}$
- $N$ is also denoted $n_{++}$ or $n_{\cdot\cdot}$

| $_\text{Row Attribute Categories}\backslash ^\text{Column Attribute Categories}$ | $1$      | $2$      | $\cdots$ | $c$      | Total    |
| ------------------------------------------------------------ | -------- | -------- | -------- | -------- | -------- |
| $1$                                                          | $n_{11}$ | $n_{12}$ | $\cdots$ | $n_{1c}$ | $n_{1+}$ |
| $2$                                                          | $n_{21}$ | $n_{22}$ | $\cdots$ | $n_{2c}$ | $n_{2+}$ |
| $\vdots$                                                     | $\vdots$ | $\vdots$ | $\ddots$ | $\vdots$ | $\vdots$ |
| $r$                                                          | $n_{r1}$ | $n_{r2}$ | $\cdots$ | $n_{rc}$ | $n_{r+}$ |
| Total                                                        | $n_{+1}$ | $n_{+2}$ | $\cdots$ | $n_{+c}$ | $N$      |

Inference and procedures are based on **conditional inference**: we treat row and column marginals *as if they were fixed*. Here "conditional" means we condition on the row/column totals. Note that sometimes one (or both) of the marginals *will* be fixed by the study design, but not always.

## Two by two table
We start with $2 \times 2$ tables for ease of discussion. Assume for now that the row variable is the `explanatory variable`, and the column variable is the `response variable`. The row variable often will have marginal totals that are fixed by design. Suppose the response variable has two levels that are "Success" and "Failure". We have

| $_\text{Explanatory} \backslash ^\text{Outcome}$ | Response 1 (S) | Response 2 (F) | Total    |
| ------------------------------------------------ | -------------- | -------------- | -------- |
| Treated                                          | $n_{11}$       | $n_{12}$       | $n_{1+}$ |
| Non-Treated                                      | $n_{21}$       | $n_{22}$       | $n_{2+}$ |
| Total                                            | $n_{+1}$       | $n_{+2}$       | $N$      |

where $n_{1+}$ and $n_{2+}$ are fixed by design, whereas $n_{+1}$ and $n_{+2}$ are not fixed.

We saw before that this design implies independent `Binomials` - $\text{Bin}(n_{1+}, P_1)$ for row 1 and $\text{Bin}(n_{2+}, P_2)$ for row 2. Our question is is there an association between the treatment and response, or in the hypothesis testing language $H_0: P_1 = P_2$, which says the probability of the successful outcome is the same for treated and non-treated.

We haven't said anything (yet) about the **common** probability of success under $H_0$, $\pi$.

In our example, we observe $n_{11}$ successes for row $1$ (treated), and $n_{21}$ successes for row $2$ (non-treated). We can use this to derive

$$
\begin{gather*}
    P_{11} = P(n_{11} \text{ successes}) = \binom{n_{1+}}{n_{11}}P_1^{n_{11}} (1 - P_1)^{n_{1+} - n_{11}} \\\\
    P_{21} = P(n_{21} \text{ successes}) = \binom{n_{2+}}{n_{21}}P_2^{n_{21}} (1 - P_2)^{n_{2+} - n_{21}} \\\\
    P_{12} = 1 - P_{11} \\\\
    P_{22} = 1 - P_{21}
\end{gather*}
$$


Also, in total we have $N$ observations with $n_{+1}$ successes:

$$
P_{+1} = \binom{N}{n_{+1}}\pi^{n_{+1}}(1-\pi)^{n_{+2}}
$$


Using properties of conditional, independent, etc. probabilities, we can show the probability of observing $n_{11}$ successes in $n_{1+}$ trails, conditional on a total of $n_{+1}$ successes in $N$ trails is

$$
P{^\ast} = \frac{P_{11}P_{21}}{P_{+1}}
$$


Under $H_0$, $\pi = P_1 = P_2$. All the probabilities in the individual Binomials end up cancelling each other out, and only the cell frequencies remain, which is the `hypergeometric distribution`:

$$
\begin{aligned}
    P^{\ast} &= \frac{(n_{11} + n_{12})!(n_{21} + n_{22})! + (n_{11} + n_{21})! + (n_{12} + n_{22})!} {(n_{11}+n_{12}+n_{21}+n_{22})!n_{11}!n_{12}!n_{21}!n_{22}!} \\\\
    &= \frac{n_{1+}!n_{2+}!n_{+1}!n_{+2}!}{N!n_{11}!n_{12}!n_{21}!n_{22}!}
\end{aligned}
$$


What might the cell counts be "expected" to look like under $H_0$? If we knew $\pi$, we could get exact expected counts for each cell using properties of the binomial. We don't know $\pi$, but we can estimate it using the sample proportion of successes: $\hat{\pi} = \frac{n_{+1}}{N}$. The expected number of $m_{11}$ (successes in row 1) can be get using $Bin(n_{1+}, \pi)$:
$$
E\left[m_{11}\right] = n_{1+}\hat{\pi} = \frac{n_{1+}n_{+1}}{N}
$$
In general,
$$
m_{ij} = \frac{\text{row i total} \times \text{col j total}}{\text{overall total}}
$$


### Odds ratio
Another useful construct is called the `odds ratio`, which is used a lot in medical settings. It can also be used as a measure of **association** between the treatment and response in a $2 \times 2$ table.

$$
\text{odds ratio} = \theta = \frac{p_1 / (1 - p_1)}{p_2 / (1 - p_2)} = \frac{\text{odds of success in row 1}}{\text{odds of success in row 2}}
$$


Under $H_0$ (no treatment effect), we saw last time that $p_1 = p_2 \Rightarrow \theta = 1$. Deviations of the sample odds ratio from 1 are indications **against** $H_0$. We can also show that

$$
\begin{gather*}
    \frac{p_1 / (1 - p_1)}{p_2 / (1 - p_2)} = \frac{m_{11}m_{22}}{m_{12}m_{21}} \\\\
    \Rightarrow \hat{\theta} = \frac{n_{11}n_{22}}{n_{12}n_{21}}
\end{gather*}
$$


### Two response variables
A second possibility for how the table could have arise: both attributes are some sort of response, e.g. both attributes could be side effects experienced when patients are given a treatment.

| $_\text{Attribute 1} \backslash ^\text{Attribute 2}$ | Level 1  | Level 2  | Total    |
| ---------------------------------------------------- | -------- | -------- | -------- |
| Level 1                                              | $n_{11}$ | $n_{12}$ | $n_{1+}$ |
| Level 2                                              | $n_{21}$ | $n_{22}$ | $n_{2+}$ |
| Total                                                | $n_{+1}$ | $n_{+2}$ | $N$      |

Only the overall sample size $N$ is fixed by design now. Each of the $N$ individuals gets allocated independently of the others to the appropriate cell $(i, j)$ according to responses on the two attributes.

Now each individual has four possible outcomes, corresponding to the different patterns of response $(y, y)$, $(y, n)$, $(n, y)$ and $(n, n)$. What this means is that all cell counts $n_{11}$, $n_{12}$, $n_{21}$, $n_{22}$ vary freely. This implies with unknown probabilities $p_{ij}$, each of the $N$ subjects is allocated independently to cell $(i, j)$, which is a `multinomial model` (multiple outcomes) with $N$ observations and probabilities $(p_{11}, p_{12}, p_{21}, p_{22})$ such that $p_{11} + p_{12} + p_{21} + p_{22} = 1$.

$p_{1+}$ and $p_{2+}$ are marginal probabilities that units fall in row $1$ or $2$, respectively; likewise, $p_{+1}$ and $p_{+2}$ are marginal probabilities that units fall in column $1$ or $2$, respectively. If the row and column attributes are independent, then $p_{ij} = p_{i+} \cdot p_{+j}$ for $i, j = 1, 2$. Since all of the probabilities here are unknown, we estimate them using the sample proportions:

$$
\hat{p}\_{i+} = \frac{n\_{i+}}{N}, \quad \hat{p}\_{+j} = \frac{n\_{+j}}{N}
$$


Under independence,
$$
\hat{p}_{ij} = \frac{n\_{i+} \cdot n\_{+j}}{N^2}
$$


By properties of the multinomial, the expected count in cell $(i, j)$ under independence is estimated to be

$$
\begin{equation} \label{eq:expected-count}
m_{ij} = N\hat{p}\_{ij} = \frac{n\_{i+} \cdot n\_{+j}}{N}
\end{equation}
$$


This is the same as in [the first model](#two-by-two-table), even though our assumptions about how the table came about were different. Also the odds ratio estimation will be the same.

### Fixing nothing
A third possibility is we collect data for $6$ months, and for each individual we classify, independently of others, to cell $(i, j)$ according to values on the two attributes. **Nothing** here is fixed, as we don't know how many individuals we're going to have!

The model here is independent `Poissons` for each of the four cells, with mean $m_{ij}$ for cell $(i, j)$. Under $H_0$ of no association between the row and column attributes, we get to the same end point as in the previous two cases.

Note that this all extends pretty easily to a general $r \times c$ table. In general, we can ignore which elements, if any, are fixed by design, and present a unified framework for the rest of the analysis.

## Unified framework for general tables
The unified framework for $r \times c$ tables focuses on level of measurement for categorical variables. We treat nominal attributes and ordinal attributes differently here.

### Nominal attributes
Three approaches are commonly used when row and column categories are both nominal: the `Fisher's exact test`, the `Pearson chi-squared test`, and the `likelihood-ratio test`. All are tests for a null of no association (or independence) between the two attributes. Under that null, all three tests have the same **large-sample** (asymptotic) distribution : $\chi^2_{(r-1)(c-1)}$. Large values are evidence against $H_0$.

Importantly, the values of the three test statistics will differ on the same sample. So will their **exact distributions**. However, unless the choice of the significance level cutoff is critical, the three seldom lead to different conclusions.

#### Fisher's exact test
We saw before for the $2 \times 2$ table in lots of detail. We also saw that for a general $r \times c$ table, "extreme" is a little murkier. We won't do it by hand. In R, the function `fisher.test()` does the job. Input either
- `x` - a two-dimensional contingency table, or
- `x` and `y` - two vectors (observations for each attribute), from which the table can be built.

The function computes **exact** p-values for the $r \times c$ table. We can also obtain a **simulated** p-value by passing the argument `simulate.p.value = T`. For $2 \times 2$ tables, it can also make inference for the odds ratio.

In practice, the Fisher's exact test is often used when asymptotic theory is inappropriate (usually when sample sizes are small), or in the case of $2 \times 2$ tables.

#### Pearson's chi-squared test
An alternative statistic for testing independence of row and column categories is the Pearson chi-squared statistic:

$$
\chi^2 = \sum_{i, j} {\frac{(m_{ij} - n_{ij})^2}{m_{ij}}}
$$


where $m_{ij}$ is the expected cell count under $H_0$ of independence, and $n_{ij}$ is the observed cell count in data. $m\_{ij}$ can be calculated using Equation $\eqref{eq:expected-count}$.

A caveat is that if there are $m_{ij} < 5$, the asymptotic $\chi^2_{(r-1)(c-1)}$ reference distribution is not a good approximation.

#### Likelihood-ratio test
An alternative statistic for test of association is the likelihood ratio:

$$
G^2 = 2 \sum_{i,j}{n_{ij} \cdot \log\left(\frac{n_{ij}}{m_{ij}}\right)}
$$


where once again $m_{ij}$ is the expected cell count under $H_0$ of independence, and $n_{ij}$ is the observed cell count in data.

On the same data, these three will obviously give different numerical results. Note that since the three tests apply to nominal categories, reordering the rows or columns doesn't affect the values of the test statistics. For ordered categories, there are more appropriate tests.

### Ordinal attributes
The following cases are considered here:
1. Nominal explanatory and ordered response.
2. Ordered explanatory and ordered response, e.g. increasing dose of a drug and varied side-effects.
3. Row and column attributes are both ordered responses, e.g. two different side-effects.

#### Nominal row and ordered column
We're going to look at the number of patients receiving each drug who experience different levels of side effects:

| $_\text{Nominal explanatory}\backslash^\text{Level (ordered response)}$ | None | Slight | Moderate | Severe | Fatal | Total |
| ------------------------------------------------------------ | ---- | ------ | -------- | ------ | ----- | ----- |
| Drug A                                                       | 23   | 8      | 9        | 3      | 2     | 45    |
| Drug B                                                       | 42   | 8      | 4        | 0      | 0     | 54    |
| Total                                                        | 65   | 16     | 13       | 3      | 2     | 99    |

> Here we're not saying anything about how this study was designed. We don't have to think about if any margins are fixed.

**Question**: is there an association between drug type (A vs. B) and level of side effect, or are these two attributes independent?

**Procedure**: This is like a massive [WMW situation]({{< ref "/courses/nonparam-stat/3-multiple-samples/3.2-two-independent-samples/index.md#test-statistic-formulations" >}}) with **lots** of ties! Both test statistic formulations can be applied.

In the Wilcoxon formulation we may use mid-ranks again, e.g. 65 subjects had no side effects, so all get the same mid-rank of 33. 16 are tied at "slight", so we set $\\{66, 67, \cdots, 81\\}$ to 73.5. We can proceed and correct for ties. 

In the Mann-Whitney formulation, there is no need to specify the mid-ranks. We just count the number of drug B patients showing the same or more severe side effects than those for each recipient of drug A, counting ties as 0.5. For instance, the 42 drug B patients with no side effects are tied with the 23 drug A patients with no side effect, $42 \times 0.5 = 21$; there are $8$ drug B patients with slight and $4$ with moderate side effects, so this group contributes $23 \cdot (42 \times 0.5 + 8 + 4) = 759$ to the test statistic. Repeat for the other responses in drug B patients.

If the nominal explanatory variable has more than two values (e.g. 3 drugs in the previous example), the obvious extension to the [Kruskal-Wallis test]({{< ref "/courses/nonparam-stat/3-multiple-samples/3.3-three-or-more-samples/index.md#the-kruskal-wallis-test" >}}) applies.

#### Ordered row and column
The row categories are ordinal explanatory, and the columns are ordered responses. We want to ask if there's an association between the two attributes. This can again be applied with the `J-T test`, with each row as an ordered sample.

**Question:** We have data on side effects experienced at increasing dose levels of a drug. Does side effects increase with dose level?

| $_\text{Dose}\backslash\text{Side effects}$ | None | Slight | Moderate | Severe |
| ------------------------------------------- | ---- | ------ | -------- | ------ |
| 100mg                                       | 50   | 0      | 1        | 0      |
| 200mg                                       | 60   | 1      | 0        | 0      |
| 300mg                                       | 40   | 1      | 1        | 0      |
| 400mg                                       | 30   | 1      | 1        | 2      |

**Procedure**: again, we score relevant ties in any column (as with the previous analysis) as $\frac{1}{2}$. Obviously the first column contributes most of the ties. We compute each pairwise `MW statistic` and add them all up. That is,

$$
\begin{aligned}
  (100, 200) &\Rightarrow 50 \cdot (\frac{60}{2} + 1) \\\\
  (100, 300) &\Rightarrow 50 \cdot (\frac{40}{2} + 2) + \frac{1}{2} \\\\
  (100, 400) &\Rightarrow 50 \cdot (\frac{30}{2} + 4) + \frac{1}{2} + 2 \\\\
  (200, 300) &\Rightarrow 60 \cdot (\frac{40}{2} + 2) + \frac{1}{2} + 1 \\\\
  (200, 400) &\Rightarrow 60 \cdot (\frac{30}{2} + 4) + \frac{1}{2} + 1 + 2 \\\\
  (300, 400) &\Rightarrow 40 \cdot (\frac{30}{2} + 4) + \frac{1}{2} + 1 + 2 + \frac{1}{2} + 2 \\\\
  \text{Total} &= 6834
\end{aligned}
$$


**Conclusion**: this has a two-sided approximate p-value of 0.4576 - no evidence of association. Even the one-sided p-value is pretty high. At first glance their results may seem counter-intuitive, but side effects are very rare so it's hard to discover much pattern.

#### The Goodman-Kruskal statistic
What if both row and column attributes are "responses"? Do high responses in one classification tend to be associated with high responses in the other (positive association) or low responses (negative association)? Or maybe there is no association between the two responses?

We can use the `J-T test` here as well, of course. The problem with the `J-T test` is that it isn't calibrated as we'd like for a measure of association. There has been a lot of work in this general class of problems - calibration for measures of association in general $r \times c$ tables (because $r$ and $c$ are relevant). Many measures have been proposed and studied.

One measure is the `Goodman-Kruskal gamma statistic` $G$. This counts concordances and discordances between row and column classifications, like Kendall's $t_k$, but with no allowance for ties. When rows and columns are ordered, for **each** count in cell $(i, j)$, there is concordance between that cell count and any cell count below and to the right:

| $_\text{Attribute A}\backslash^\text{Attribute B}$ | $\cdots$           | $B_j$        | $\cdots$        |
| -------------------------------------------------- | ------------------ | ------------ | --------------- |
| $\vdots$                                           | $\vdots$           | $\vdots$     | $\vdots$        |
| $A_i$                                              | $\leftarrow$       | $(i, j)$     | $\rightarrow$   |
| $\vdots$                                           | ordering disagrees | $\downarrow$ | ordering agrees |

Cell $(i, j)$ has count $n_{ij}$. Let $N_{ij}^+$ denote the sum of all counts below and to the right of cell $(i, j)$, not including the ones in the same row/column. The total number of condordances $C$ is given by

$$
C = \sum_{i, j}{n_{ij} \cdot N_{ij}^+} \qquad 1 \leq i \leq r - 1, \quad 1 \leq j \leq c - 1
$$

Similarly, let $N_{ij}^-$ denote the sum of all counts in cells to the left and below cell $(i, j)$, which is ordering on the two attributes disagrees. The total number of discordances $D$ is given by
$$
D = \sum_{i, j}{n_{ij} \cdot N_{ij}^-} \qquad 1 \leq i \leq r - 1, \quad 2 \leq j \leq c
$$


We define our test statistic $G$:

$$
G = \frac{C-D}{C+D}
$$


which is calibrated to be between $-1 (C = 0)$ and $+1(D = 0)$ so that we can easily interpret.

## Testing goodness of fit
> Here we'll only talk about the $\chi^2$ test, where we have an asymptotic $\chi^2$ distribution as our reference. It's widely used as a goodness-of-fit test of data to any discrete distribution.

Instead of a model of independence or lack of association as before, we can consider goodness of fit to a hypothesized discrete distribution, which may be binomial, Poisson, uniform or some other discrete distribution. We compute **expected cell counts** under the hypothesized model, and compare to the **observed** counts in the data.

Sometimes, the parameter(s) of those distributions will not be known. In this case, we need to **estimate** the unknown parameters from the data.

### All parameters are known
A computer program is supposed to generate random digits from 0 to 9. If it is doing so, we'll get digits that look like i.i.d. observations on the values 0 to 9, each with probability 0.1.

We want to test $H_0$: the numbers appear to be random digits, vs. $H_1$: some digits are more likely than others.

Generate some number, suppose 300, of digits from the program, and compare to the discrete uniform on 0 to 9 under $H_0$. Under $H_0$, the **expected** number of each number is $30$. The observed results in our sample:

| Number | Observation |
| ------ | ----------- |
| 0      | 22          |
| 1      | 28          |
| 2      | 41          |
| 3      | 35          |
| 4      | 19          |
| 5      | 25          |
| 6      | 25          |
| 7      | 40          |
| 8      | 30          |
| 9      | 35          |

$$
\begin{aligned}
    \chi^2 &= \sum_{i=0}^9{\frac{(O_i - E_i)^2}{E_i}} \\\\
    &= \frac{(22-30)^2}{30} + \frac{(28-30)^2}{30} + \cdots + \frac{(35-30)^2}{30} \\\\
    &= 17
\end{aligned}
$$



Compare this to $\chi^2_{(9)}$ ($r$ categories has $r-1$ degrees of freedom for $\chi^2$).

In R, use the `chisq.test()` function. In this case:

```r
x <- c(22, 28, 41, 35, 19, 25, 25, 40, 30, 35)
chisq.test(x)
```

The default is to test for **uniform** if no other parameters are specified. This gives a p-value of approximately 0.049. For a different set of specified probabilities, e.g. as given by binomial or Poisson, we need to pass another vector of the same length as the data vector specifying these probabilities. Alternatively, we can supply the vector of expected counts under $H_0$ and specify `rescale.p = T`.

### Test with estimated parameters
If we have to estimate parameters, we lose a degree of freedom for each one. In this situation, R will compute the test statistic but won't use the correct degrees of freedom. Use `pchisq()` with the correct degrees of freedom in this case.

Suppose we want test the goodness of fit to a binomial, but the probability of success $P$ is unknown and we have to estimate it.

We have data on the first 18 major league baseball players to have 45 times at bat in 1970. The number of hits they got in their 45 times at bat are given as follows:

$$
\begin{gather*}
  18 && 17 && 16 && 15 && 14 && 14 && 13 && 12 && 11 \\\\
  11 && 10 && 10 && 10 && 10 && 10 && 9 && 8 && 7
\end{gather*}
$$


We will test the null hypothesis that these data follow a $Bin(45, p)$. But first we need to estimate $p$ for each time at bat. We may use the overall relative frequency across all players of getting a hit: 

\\[\hat{p} = \text{proportion of hits out of the total times at bat} = 0.2654\\]

With $n=45$ and $\hat{p} = 0.2654$, we then compute the probability of $0, 1, \cdots, 45$ hits from $Bin(45, 0.2654)$ and get expected counts. Theoretically we should get all expected counts from 0 to 45, but there are many small values which will result in $0$ expected counts.

```r
head(18 * dbinom(0:45, size = 45, prob = 0.2654))
# [1] 1.688745e-05 2.745533e-04 2.182224e-03
# [4] 1.130047e-02 4.286825e-02 1.269988e-01
```

Here we multiply by 18 because there are 18 players. We get **many** small probabilities, which means we'd get many expected counts of $0$, and this is a problem for $\chi^2$ tests! In practice, we should combine some values to get more reasonable expected cell counts, e.g. combine cells with expected counts less than $0.5$, so the first $8$ values, and all values $18$ and above are combined.

```r
# 7 hits or below
sum(18 * dbinom(0:45, size = 45, prob = 0.2654)[1:8])
# [1] 1.105234

18 * dbinom(0:45, size = 45, prob = 0.2654)[9:18]
# [1] 1.0566192 1.5693786 2.0411749 2.3464190 2.4018906
# [6] 2.2027936 1.8190547 1.3582077 0.9200627 0.5670437

# 18 hits or above
sum(18 * dbinom(0:45, size = 45, prob = 0.2654)[19:46])
# [1] 0.6121209
```

Now we can build the chi-square goodness-of-fit test:

```r
obs <- c(1, 1, 1, 5, 2, 1, 1, 2, 1, 1, 1, 1)
expected.values <- c(
  round(sum(18 * dbinom(0:45, size = 45, prob = 0.2654)[1:8]), 2),
  round(18 * dbinom(0:45, size = 45, prob = 0.2654)[9:18], 2),
  round(sum(18 * dbinom(0:45, size = 45, prob = 0.2654)[19:46]), 2)
)
chisq.test(obs, p = expected.values, rescale.p = T)

# 	Chi-squared test for given probabilities
#
# data:  obs
# X-squared = 6.737, df = 11, p-value = 0.82
#
# Warning message:
# In chisq.test(obs, p = expected.values, rescale.p = T) :
#   Chi-squared approximation may be incorrect
```

The value of the test statistic is correct, but we should have 10 degrees of freedom instead of 11, because estimating $\hat{p}$ from the data costs us 1 degree of freedom!

```r
pchisq(6.737, df = 10, lower.tail = F)
# [1] 0.7500185
```

Here we specify `lower.tail = F` because we want the probability to the right of our observed test statistic value, and by default the probability to the left is calculated.

---

This is the last part of the "classical" nonparametric statistics. Next, we'll be focusing on topics in modern nonparametric statistics, which is also the finale of our nonparametric methods discussions.