---
title: "Bayesian Linear Regression"
date: 2021-04-05T10:20:00-05:00
summary: "The main difference with traditional approaches is in the specification of prior distributions for the regression parameters, which relate covariates to a continuous response variable. However, the Bayesian approach also provides a fairly intuitive way to add random effects (such as a random intercept or random slope), which results in what is traditionally known as a linear mixed model." # appears in list of posts
categories: ["Bayesian Statistics"] # main category; shown in post metadata
tags: ["Statistics", "Bayesian Statistics", "Visualization", "Estimation"] # list of related tags

slug: "bayesian-stat-bayesian-linear-regression"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 100 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Curvy road in Sedona, Arizona." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "DlnK1KOREds" # Unsplash ID of the picture
---



Unlike many of the previous lectures, a lot of this might not feel very different from traditional approaches because we are placing many uninformative priors. However, the Bayesian approach does have some benefits, which we will touch upon as we go through 3 different regression models on the following data set.


## London schools example

The `student.txt` data set from Goldstein et al.[^goldstein] contains externally standardized results of 1978 students' performance on the General Certificate of Secondary Examination (GCSE), an exam taken by British students finishing Year 11 (about age 15). These students were drawn from 38 schools within inner London. In addition, the following information is included in the data:

- Standardized scores on the London Reading Test (LRT), an examination required for all British students at age 11.
- Students were rated into one of three groups based on a Verbal Reasoning (VR) test taken at age 11 -- `\(\text{VR}=1\)` being the best score and `\(\text{VR}=3\)` being the worst score.
- Gender at birth coded 0 if male and 1 if female.

[^goldstein]: Goldstein, H., Rasbash, J., Yang, M., Woodhouse, G., Pan, H., Nuttall, D., & Thomas, S. (1993). A Multilevel Analysis of School Examination Results. *Oxford Review of Education, 19*(4), 425-433. Retrieved April 7, 2021, from http://www.jstor.org/stable/1050563

Our goal is to model the relationship between GCSE performance and these covariates. With a relatively small number of covariates, we can examine their pairwise relationships with the response variable. We begin by ignoring which school each student attends for now, and plot GCSE against LRT.


```r
library(tidyverse)
library(ggpubr)
library(rjags)
library(patchwork)
```


```r
dat <- data.table::fread("student.txt", sep = "\t") %>%
  select(Y, LRT, VR1, VR2, Gender, School)

ggscatter(dat, x = "LRT", y = "Y", ylab = "GCSE", size = 1)
```

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/gcse-vs-lrt-1.svg" alt="Scatter plot of GCSE vs. LRT." width="576" />
<p class="caption">Figure 1: Scatter plot of GCSE vs. LRT.</p>
</div>

Both variables are centered at zero. It seems in general a higher LRT score corresponds to a higher GCSE score.

Next we have the gender at birth variables. The GCSE scores by gender are summarized in the following table.


|Gender | Sample size|    Mean| Standard deviation|
|:------|-----------:|-------:|------------------:|
|Female |        1209|  0.0926|             0.9700|
|Male   |         769| -0.1453|             1.0351|

There are fewer males than females in the sample, and it seems like females outperform the average and the males tend to perform below par.

The final variable is VR, and we observe similar trends where `\(\text{VR}=1\)` is the best-performing group, whereas `\(\text{VR}=3\)` performed the worst.


|VR | Sample size|    Mean| Standard deviation|
|:--|-----------:|-------:|------------------:|
|1  |         521|  0.8101|             0.8716|
|2  |        1160| -0.1310|             0.8426|
|3  |         297| -0.9083|             0.7385|


## Model 1 - multiple linear regression

Same as always, we need to first specify our **sampling model**. This part of the model relates `\(Y_{ij}\)` to the explanatory variables. Let `\(Y_{ij}\)` be the values of the GCSE for student `\(i\)` within school `\(j\)`:

<div>
$$
Y_{ij} \mid \boldsymbol{\beta}, \sigma^2 \overset{\text{indep.}}{\sim} \mathcal{N}(\mu_{ij}, \sigma^2), \quad i = 1, \cdots, n_j,\, j = 1, \cdots, J
$$
</div>

where `\(\boldsymbol{\beta}\)` is a collection of five regression coefficients `\((\beta_1, \cdots, \beta_5)\)`, and the "regression function" is:

<div>
$$
\mu_{ij} = \beta_1 + \beta_2(\text{LRT})_{ij} + \beta_3 (\text{VR1})_{ij} + \beta_4 (\text{VR2})_{ij} + \beta_5 (\text{Gender})_{ij}
$$
</div>

Assuming conditional independence across the samples, the probability density function is:

<div>
$$
p(y \mid \boldsymbol{\beta}, \sigma^2) = \prod_{j=1}^J \prod_{i=1}^{n_j} p(y_{ij} \mid \boldsymbol{\beta}, \sigma^2)
$$
</div>

Next up is the **prior model**. The unknown values in the sampling model are `\(\mu_{ij}\)` and `\(\sigma^2\)`, but note that the mean function `\(\mu_{ij}\)` is determined by the regression coefficients `\(\beta\)`'s, so we need to specify priors for six parameters.

Suppose we don't have any information about the parameters, then:

<div>
$$
\left.\begin{aligned}
  \sigma^2 \sim \text{IG}(0.01, 0.01)& \\
  \beta_1, \cdots, \beta_5 \overset{iid}{\sim} N(0, 10^{10})&
\end{aligned}\right\} \text{independent}
$$
</div>

We've seen the inverse-gamma prior multiple times for the variance term. For the normal priors, they are proper priors but are virtually "flat", so it's almost like putting a uniform prior. The probability density function is:

<div>
$$
p(\boldsymbol{\beta}, \sigma^2) = \left( \prod_{i=1}^5 p(\beta_i) \right) p(\sigma^2)
$$
</div>

We cannot write down the full **posterior distribution** in closed-form. However, for our prior specification, it can be shown that the full conditional distributions do exist in closed-form. The ramifications of that is a *Gibbs sampler* can be used to generate approximate samples from the posterior[^jags]. 
[^jags]: JAGS will recognize this automatically. In Bayesian sampling, a Gibbs sampler is usually preferable to the Metropolis-Hastings algorithm, because the latter need to compute the acceptance ratio and is thus slower.


### Posterior samples via JAGS

Now we can draw posterior samples using JAGS. Again, the model specified below is without the effect of schools.


```r
## (1) Create data list
n <- nrow(dat)      # sample size
p <- ncol(dat) - 1  # number of regression coefficients (-school-Y+intercept)

dataList1 <- list(
  "n" = n,
  "p" = p,
  "Y" = dat$Y,
  "LRT" = dat$LRT,
  "VR1" = dat$VR1,
  "VR2" = dat$VR2,
  "Gender" = dat$Gender
)

## (2) Specify list of parameter(s) to be monitored
parameters1 <- c("beta", "sig2")

## (3) Specify initial values for parameter(s) in Metropolis-Hastings algorithm
initsValues1 <- list(
  "beta" = rep(0, p),
  "tau2" = 1
)

## (4) Specify parameters for running Metropolis-Hastings algorithm
adaptSteps <- 1000     # number of steps to "tune" the samplers
burnInSteps <- 2000    # number of steps to "burn-in" the samplers
nChains <- 3            # number of chains to run
numSavedSteps <- 5000  # total number of steps in chains to save
thinSteps <- 10         # number of steps to "thin" (1 = keep every step)
nIter <- ceiling((numSavedSteps*thinSteps)/nChains) 	# steps per chain

## (5) Create, initialize, and adapt the model
# This will require you to create a separate .txt file which specifies
# the model
model1 <- textConnection("
model
{
  # Likelihood - in JAGS, normal distribution is parameterized by
  # mean theta and precision = tau2 = 1/sig2
  for (i in 1:n) {
  	Y[i] ~ dnorm(mu[i], tau2)
    mu[i] = beta[1]+beta[2]*LRT[i]+beta[3]*VR1[i]+beta[4]*VR2[i]+beta[5]*Gender[i]
  }

  # Priors
  for (i in 1:p) {
    beta[i] ~ dnorm(0, 1e-10)
  }

  tau2 ~ dgamma(0.01, 0.01)

  # Need to have model calculate sig2 = 1/tau2
  sig2 <- 1/tau2
}
")
jagsModel1 <- jags.model(model1, 
                         data = dataList1, 
                         inits = initsValues1, 
                         n.chains = nChains, 
                         n.adapt = adaptSteps)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 1978
##    Unobserved stochastic nodes: 6
##    Total graph size: 10236
## 
## Initializing model
```

Now we actually draw the samples (takes quite a while):


```r
## (6) Burn-in the algorithm
if (burnInSteps > 0) {
  update(jagsModel1, n.iter = burnInSteps)
}

## (7) Run MCMC algorithm
codaSamples1 <- coda.samples(jagsModel1, 
                             variable.names = parameters1, 
                             n.iter = nIter, 
                             thin = thinSteps)
close(model1)
# summary(codaSamples1)
```

The posterior summaries of `\(\beta\)` and `\(\sigma^2\)` are given below. First thing we might notice is VR1 (`\(\beta_3\)`) appears to have the strongest positive connection with GCSE scores. VR2 also has a slight positive effect.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/model1-params-boxplot-1.svg" alt="Posterior samples of betas." width="672" />
<p class="caption">Figure 2: Posterior samples of betas.</p>
</div>

The gender (`\(\beta_5\)`) boxplot is primarily above zero, which indicates females generally have higher GCSE scores. LRT (`\(\beta_2\)`) also appears to have a very small but positive relationship with GCSE. Recall that LRT spanned (-30, 30), and we should consider standardizing the variable to get a better sense of the relationship between LRT and GCSE.

Finally, the intercept indicating that the baseline level (LRT=0, Gender=0, VR3), is negative.


### Notes

As with traditional linear regression methods, we should consider the following items.

1. Is the goal of the model explanatory or prediction? If explanatory, are there any predictor variables which are redundant or provide minimal information in explaining changes in the response? 
2. Consider transformations if assumptions do not appear to be satisfied.
    - Standardization of the predictor variables.
    - If GCSE wasn't normally distributed, maybe try log(GCSE).
    - Other transformations of the predictors if we observe improved correlation with the response in the EDA.

3. Including higher-order interaction terms. For example,

    <div>
    $$
    \mu_{ij} = \beta_1 + \beta_2 x_1 + \beta_3 x_2 + \beta_4 x_1 x_2
    $$
    </div>

There's of course other things to consider (e.g. regularization), and we'll talk about some other types of priors in the next lecture.


## Model 2 - adding a school effect

Our first linear regression model ignores which school each student attended. Should we account for a school effect in our model? We might be interested in which schools are on average performing better when the other variables are factored out.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/gsce-school-boxplot-1.svg" alt="Boxplot of GCSE vs. School." width="960" />
<p class="caption">Figure 3: Boxplot of GCSE vs. School.</p>
</div>

The boxplot shows some variability across the schools. Some are almost centered around 1 indicating high-performance. Note that there's also varying sample sizes, and some schools only have one observation.

So how could Model 1 be modified to account for school-specific average GCSE performance? We may specify the **sampling model** as:

<div>
$$
Y_{ij} \mid \boldsymbol{\alpha}, \boldsymbol{\beta}, \sigma^2 \overset{\text{indep.}}{\sim} \mathcal{N}(\mu_{ij}, \sigma^2)
$$
</div>

where `\(\boldsymbol{\alpha} = (\alpha_1, \cdots, \alpha_J)\)`, `\(\boldsymbol{\beta} = (\beta_1, \cdots, \beta_5)\)`, and

<div>
$$
\mu_{ij} = \alpha_j + \beta_1 + \beta_2(\text{LRT})_{ij} + \beta_3 (\text{VR1})_{ij} + \beta_4 (\text{VR2})_{ij} + \beta_5 (\text{Gender})_{ij}.
$$
</div>

This is very similar to Model 1, except in the regression function we now have an `\(\alpha_j\)` term that models school-specific effect for school `\(j\)`.

Under this new model, we assume that the predictors have the same relationship with the response variable *regardless* of the school. This is called a `random intercept model`, where the slopes are the same but the intercepts differ.

Next we specify **priors** for all our parameters. For the `\(\beta\)`'s we would use exactly the same prior as in Model 1:

<div>
$$
\beta_1, \cdots, \beta_5 \overset{iid}{\sim} \mathcal{N}(0, 10^{10})
$$
</div>

For the `\(\alpha\)`'s, since the schools were randomly selected from inner London, we would expect them to share some similarities after accounting for all other variables. So following the setup in the hierarchical model, the prior for the school effects could be:

<div>
$$
\alpha_1, \cdots, \alpha_J \mid \sigma_\alpha^2 \sim \mathcal{N}(0, \sigma_\alpha^2)
$$
</div>

where `\(\sigma_\alpha^2\)` models the variability between schools. Assuming the `\(\alpha\)`'s are dependent encourages shrinkage of school effects, and improves estimation by borrowing information from other schools.

We also have priors for the variance parameters:

`\begin{gathered}
  \sigma_\alpha^2 \sim \text{IG}(0.01, 0.01) \\
  \sigma^2 \sim \text{IG}(0.01, 0.01)
\end{gathered}`

We assume all of these parameters are independent. The final Model 2 is a hierarchical linear model with school-specific random intercepts, i.e. treating the intercept as a `random effect`.

### Posterior samples via JAGS

Using the code below, we can obtain posterior samples of the parameters.


```r
n.schools <- length(unique(dat$School))  # number of schools
dataList2 <- list(
  "n" = n,
  "p" = p,
  "n.schools" = n.schools,
  "Y" = dat$Y,
  "LRT" = dat$LRT,
  "VR1" = dat$VR1,
  "VR2" = dat$VR2,
  "Gender" = dat$Gender,
  "School" = dat$School
)

parameters2 <- c("alpha", "beta", "sig2", "sig2.alpha", "icc")

initsValues2 <- list(
  "alpha" = rep(0, n.schools),
  "beta" = rep(0, p),
  "tau2" = 1,
  "tau2.alpha" = 1
)

model2 <- textConnection("
model {
  # Likelihood - in JAGS, normal distribution is parameterized by
  # mean theta and precision = tau2 = 1/sig2
  for (i in 1:n) {
  	Y[i] ~ dnorm(mu[i],tau2)
    mu[i] = alpha[School[i]]+beta[1]+beta[2]*LRT[i]+beta[3]*VR1[i]+beta[4]*VR2[i]+beta[5]*Gender[i]
  }

  # Priors
  for (i in 1:p) {
    beta[i] ~ dnorm(0,1e-10)
  }

  for (j in 1:n.schools) {
    alpha[j] ~ dnorm(0,tau2.alpha)
  }

  tau2 ~ dgamma(0.01,0.01)
  tau2.alpha ~ dgamma(0.01,0.01)

  # Need to have model calculate variances
  sig2 = 1/tau2
  sig2.alpha = 1/tau2.alpha

  # Perhaps we also want the intraclass correlation
  icc = sig2.alpha/(sig2+sig2.alpha)
}
")
jagsModel2 <- jags.model(model2, 
                         data = dataList2, 
                         inits = initsValues2, 
                         n.chains = nChains, 
                         n.adapt = adaptSteps)
```

```
## Compiling model graph
##    Resolving undeclared variables
##    Allocating nodes
## Graph information:
##    Observed stochastic nodes: 1978
##    Unobserved stochastic nodes: 45
##    Total graph size: 13373
## 
## Initializing model
```

```r
if (burnInSteps > 0) {
  update(jagsModel2, n.iter = burnInSteps)
}

codaSamples2 <- coda.samples(jagsModel2, 
                             variable.names = parameters2, 
                             n.iter = nIter,
                             thin = thinSteps)
close(model2)
```

The summaries for the `\(\beta\)`'s are very similar to the ones we got in Model 1. This is because Model 1 "averaged" all the school effects. 

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/model2-beta-1.svg" alt="Posterior samples of betas." width="672" />
<p class="caption">Figure 4: Posterior samples of betas.</p>
</div>

What's more interesting is the school-specific intercepts. Any `\(\alpha\)` value that's above zero indicates a school that outperformed the baseline. School 9 appears to be the "best" school in terms of the median, and the worst-performing school could be school 17. For schools 31-38 with smaller sample sizes, their school effects are pushed towards zero and the posterior samples have a lot of variation.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/model2-alpha-1.svg" alt="Posterior samples of alphas." width="960" />
<p class="caption">Figure 5: Posterior samples of alphas.</p>
</div>

Figure <a href="#fig:model2-alpha">5</a> can be difficult to read with this many schools. An alternative is to produce an *overall ranking of schools*[^school-ranks] by averaging over the sample-to-sample rankings of school-specific effect parameters `\(\alpha_1, \cdots, \alpha_{38}\)`:

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/school-overall-ranking-1.svg" alt="School rankings where points are posterior means and lines are 95% CIs." width="768" />
<p class="caption">Figure 6: School rankings where points are posterior means and lines are 95% CIs.</p>
</div>

[^school-ranks]: R code for generating the school rankings figure:

    
    ```r
    # Used rank of mean instead of mean of ranks.
    # With the large sample size results are the same.
    mcmc2 %>%
      select(starts_with("alpha")) %>%
      pivot_longer(everything(), names_to = "School", values_to = "Value") %>%
      mutate(School = str_replace(School, "^alpha\\[(\\d+)\\]$", "\\1")) %>%
      group_by(School) %>%
      summarise(
        post_mean = mean(Value),
        lowerCI = quantile(Value, c(0.025,0.975))[1],
        upperCI = quantile(Value, c(0.025,0.975))[2]
      ) %>%
      ungroup() %>%
      arrange(post_mean) %>%
      mutate(School = fct_inorder(School)) %>%
      ggscatter(x = "post_mean", y = "School", xlab = "") +
      geom_segment(aes(x = lowerCI, xend = upperCI, y = School, yend = School)) +
      xlim(-1, 1)
    ```
    
Finally, we can compute the posterior density for the `intraclass correlation` defined by the quantity

<div>
$$
\frac{\sigma_\alpha^2}{\sigma^2 + \sigma_\alpha^2}
$$
</div>

to estimate the proportion of total variation due to school-specific effects. This is of interest because the schools are a random sample of a population of schools. If the ICC is high, then it means the grouping explains a lot of the total variation. 

<img src="{{< blogdown/postref >}}index_files/figure-html/intraclass-var-1.svg" width="480" />

The ICC is centered around 0.08, which mean about 8% of the total variation in GCSE scores is attributed to school effects, and the remaining 92% is just associated with variation from student to student. We care about quantities like this more in random effect models (vs. fixed effect models) because we can generalize this and infer about all schools instead of just those 38 schools.


## Model 3 - adding school-level covariates

The final question we may want to consider is what factors about *the schools themselves* drive GCSE performance? For example, what factors attributed to school 9 influenced high GCSE performance?

A separate `school.txt` data set contains information about each of the 38 schools. Each school belongs to one of four *denominations* -- public (baseline category), church of England (CE), Roman Catholic (RC), or other. The *School gender* is either mixed gender (baseline category), girls only, or boys only.


| School| CE| RC| Other| Girls| Boys|
|------:|--:|--:|-----:|-----:|----:|
|      1|  0|  0|     0|     0|    0|
|      2|  1|  0|     0|     0|    0|
|      3|  0|  0|     0|     0|    0|
|      4|  0|  0|     0|     0|    0|
|      5|  0|  0|     0|     0|    1|
|      6|  0|  0|     0|     0|    0|
|      7|  0|  0|     0|     0|    1|
|      8|  0|  1|     0|     0|    0|
|      9|  0|  1|     0|     0|    0|
|     10|  0|  1|     0|     0|    1|
|     11|  0|  0|     1|     0|    0|
|     12|  0|  0|     0|     0|    0|
|     13|  0|  0|     0|     1|    0|
|     14|  0|  0|     0|     0|    1|
|     15|  0|  0|     0|     0|    0|
|     16|  0|  0|     0|     1|    0|
|     17|  1|  0|     0|     0|    0|
|     18|  0|  0|     0|     0|    1|
|     19|  0|  1|     0|     0|    1|
|     20|  0|  0|     1|     0|    1|
|     21|  0|  0|     0|     0|    0|
|     22|  0|  1|     0|     0|    0|
|     23|  0|  0|     0|     0|    0|
|     24|  0|  0|     0|     0|    1|
|     25|  0|  0|     0|     1|    0|
|     26|  0|  0|     0|     1|    0|
|     27|  0|  0|     0|     0|    0|
|     28|  0|  1|     0|     0|    1|
|     29|  0|  1|     0|     1|    0|
|     30|  0|  0|     0|     0|    1|
|     31|  1|  0|     0|     0|    0|
|     32|  0|  0|     0|     0|    0|
|     33|  0|  0|     0|     0|    0|
|     34|  0|  0|     0|     0|    1|
|     35|  0|  0|     0|     0|    1|
|     36|  0|  1|     0|     0|    0|
|     37|  0|  1|     0|     0|    0|
|     38|  0|  0|     1|     0|    0|

This additional information might give us insight to help explain some of the variation in the school-specific baseline performance. Comparing GCSE scores by school denomination and school gender[^gcse-covariates], we can see public schools are somewhat average, and RC schools seem to perform best; there's not a large difference between different school gender types.

<img src="{{< blogdown/postref >}}index_files/figure-html/gcse-vs-denomination-school-gender-1.svg" width="960" />

[^gcse-covariates]: The R code for generating the boxplots comparing GCSE scores by school denomination and school gender:

    
    ```r
    dat_p <- dat %>%
      select(GCSE = Y, School) %>%
      inner_join(schools, by = "School")
    
    p1 <- dat_p %>%
      mutate(Denomination = case_when(
        CE == 1 ~ "CE",
        RC == 1 ~ "RC",
        Other == 1 ~ "Other",
        T ~ "Public"
      )) %>%
      ggboxplot(x = "Denomination", y = "GCSE")
    
    p2 <- dat_p %>%
      mutate(`School gender` = case_when(
        Girls == 1 ~ "Girls",
        Boys == 1 ~ "Boys",
        T ~ "Mixed"
      )) %>%
      ggboxplot(x = "School gender", y = "GCSE")
    
    p1 | p2
    ```

To incorporate this into our model, we can alter our school-specific level of the previous hierarchical linear regression model (Model 2). The student-level of our **sampling model** stays the same:

<div>
$$
Y_{ij} \mid \boldsymbol{\alpha}, \boldsymbol{\beta}^{(\text{st})}, \sigma^2 \overset{\text{indep.}}{\sim} \mathcal{N}(\mu_{ij}^{(\text{st})}, \sigma_\alpha^2)
$$
</div>

and the regression function is the same as in Model 2, only with additional (st) superscripts indicating student-specific effects. For the school-level sampling model:

<div>
$$
\alpha_j \mid \boldsymbol{\beta}^{(\text{sc})}, \sigma_\alpha^2 \overset{\text{indep.}}{\sim} \mathcal{N}(\mu_j^{(\text{sc})}, \sigma_\alpha^2)
$$
</div>

The difference here is instead of a normal distribution centered at zero, we have a `\(\mu_j^{(\text{sc})}\)` term which we will regress school effects on the school-specific variables:

<div>
$$
\mu_j^{(\text{sc})} = \beta_1^{(\text{sc})} (\text{CE})_j + \beta_2^{(\text{sc})} (\text{RC})_j + \beta_3^{(\text{sc})} (\text{Other})_j + \beta_4^{(\text{sc})} (\text{Girls})_j + \beta_5^{(\text{sc})} (\text{Boys})_j
$$
</div>

Note that there's **no intercept term** because we want the mean of the normal distribution to be zero.

With this, our **prior model** is:

<div>
$$
\left.\begin{align}
  \beta_1^{(\text{st})}, \cdots, \beta_5^{(\text{st})} \overset{iid}{\sim} \mathcal{N}(0, 10^{10}) \\
  \beta_1^{(\text{sc})}, \cdots, \beta_5^{(\text{sc})} \overset{iid}{\sim} \mathcal{N}(0, 10^{10}) \\
  \sigma^2 \sim \text{IG}(0.01, 0.01) \\
  \sigma_\alpha^2 \sim \text{IG}(0.01, 0.01)
\end{align}\right\} \text{independent}
$$
</div>


### Posterior samples



[^model3]: R code for specifying Model 3:

    
    ```r
    p.schools <- ncol(schools) - 1
    
    dataList3 <- list(
      "n" = n,
      "p" = p,
      "n.schools" = n.schools,
      "Y" = dat$Y,
      "LRT" = dat$LRT,
      "VR1" = dat$VR1,
      "VR2" = dat$VR2,
      "Gender" = dat$Gender,
      "School" = dat$School,
      "p.schools" = p.schools,
      "CE" = schools$CE,
      "RC" = schools$RC,
      "Other" = schools$Other,
      "Girls" = schools$Girls,
      "Boys" = schools$Boys
    )
    
    parameters3 <- c("alpha", "beta_st", "beta_sc", "sig2", "sig2.alpha")
    
    initsValues3 <- list(
      "alpha" = rep(0, n.schools),
      "beta_st" = rep(0, p),
      "beta_sc" = rep(0, p.schools),
      "tau2" = 1,
      "tau2.alpha" = 1
    )
    
    model3 <- textConnection("
    model {
      # Likelihood - in JAGS, normal distribution is parameterized by
      # mean theta and precision = tau2 = 1/sig2
      for (i in 1:n) {
      	Y[i] ~ dnorm(mu_st[i], tau2)
        mu_st[i] = alpha[School[i]]+beta_st[1]+beta_st[2]*LRT[i]+beta_st[3]*VR1[i]+beta_st[4]*VR2[i]+beta_st[5]*Gender[i]
      }
    
      for(j in 1:n.schools){
        alpha[j] ~ dnorm(mu_sc[j],tau2.alpha)
        mu_sc[j] = beta_sc[1]*CE[j]+beta_sc[2]*RC[j]+beta_sc[3]*Other[j]+beta_sc[4]*Girls[j]+beta_sc[5]*Boys[j]
      }
    
      # Priors
      for (i in 1:p) {
        beta_st[i] ~ dnorm(0, 1e-10)
      }
    
      for(j in 1:p.schools){
        beta_sc[j] ~ dnorm(0, 1e-10)
      }
    
      tau2 ~ dgamma(0.01, 0.01)
      tau2.alpha ~ dgamma(0.01, 0.01)
    
      # Need to have model calculate variances
      sig2 = 1/tau2
      sig2.alpha = 1/tau2.alpha
    }
    ")
    jagsModel3 <- jags.model(model3, 
                             data = dataList3, 
                             inits = initsValues3, 
                             n.chains = nChains, 
                             n.adapt = adaptSteps)
    close(model3)
    
    if (burnInSteps > 0) {
      update(jagsModel3, n.iter = burnInSteps)
    }
    
    codaSamples3 <- coda.samples(jagsModel3, 
                                 variable.names = parameters3, 
                                 n.iter = nIter, 
                                 thin = thinSteps)
    ```
    
With the posterior samples obtained via JAGS[^model3], we may first examine the student parameters.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/model3-beta-student-1.svg" alt="Posterior samples of beta_st." width="672" />
<p class="caption">Figure 7: Posterior samples of beta_st.</p>
</div>

Once again we get the same relationships with models 1 and 2, where VR1 is the best predictor of GCSE scores. The benefit is we've now accounted for all the other variables.

<div class="figure">
<img src="{{< blogdown/postref >}}index_files/figure-html/model3-beta-school-1.svg" alt="Posterior samples of beta_sc." width="672" />
<p class="caption">Figure 8: Posterior samples of beta_sc.</p>
</div>

The posterior summaries of the school-specific parameters `\(\beta^{(\text{sc})}\)` show that RC schools tend to have higher GCSE performance. Girls-only and boys-only schools also tend to have higher GCSE performance when all other variables are accounted for, although the difference wasn't so significant in Figure <a href="#fig:gcse-vs-denomination-school-gender"><strong>??</strong></a>. 

## Final notes

We discussed three regression models in this lecture, starting from the traditional multiple regression model (Model 1). Combining linear regression with hierarchical modeling yields Bayesian formulations of linear mixed models (Models 2 & 3). We focused on models with random intercepts, but we could have also considered models with *random slopes* or *random intercepts and slopes*.

If we have many predictors in JAGS, it will be much more convenient to use the `inprod` function when specifying the model, i.e.:

```r
Y[i] ~ dnorm(mu[i], tau2)
mu[i] <- alpha + inprod(X[i, ], beta)
```

where `X[i, ]` is the design matrix.

In the next lecture, we will talk about model selection and account for non-normal response data.
