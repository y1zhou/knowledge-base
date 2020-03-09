---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Basic Concepts"
categories:
  - Nonmarametric Methods
  - Statistics
tags:
  - Nonmarametric Methods
  - Statistics
summary: A brief introduction to what we're going to discuss in later chapters.
date: 2019-01-25T22:50:34-05:00
lastmod: 2019-01-25T22:50:34-05:00
draft: false  # Is this a draft? true/false
toc: false # Show table of contents? true/false
type: docs  # Do not modify.

featured: true
weight: 1

# Add menu entry to sidebar.
# - Substitute `example` with the name of your course/documentation folder.
# - name: Declare this menu item as a parent with ID `name`.
# - parent: Reference a parent ID if this page is a child.
# - weight: Position of link in menu.
menu:
  nonparam-stat:
    name: Basic Concepts
    parent: Introduction
    weight: 1
---

We want to move away from "standard" or "typical" approaches to statistical inference, where we assume that our data are drawn from some distributional family, e.g. the standard setup in which $X_1, X_2, ..., X_n \sim N(\mu, \sigma^2)$. Here $N(\mu, \sigma^2)$ is a normal distributional family. Similarly we could have $Pois(\lambda)$ for a Poisson distribution. In these cases, we're making assumptions about the underlying distribution. These assumptions may (or may not) be realistic or valid. In any case, they are restrictive.

Nonparametric statistical methods (sometimes called "distribution-free" methods) aim to {{< hl >}}relax these assumptions about distributional forms{{</hl>}}. They will be more **general** and more **robust** [^1], but we sacrifice **power** (not always) if the data truly come from a particular family, such as Normal, for which optimal tests (such as `z-test` or `t-test`) exist.

The term `nonparametric method` is also used in a variety of ways, which we want to examine:

- Classical approaches, e.g. based on `ranks`
- Computational approaches, e.g. `bootstrap`
- Modern regression (and other) approaches, e.g. `smoothing`

The big question is if we don't assume a distributional family, how can we proceed to do inference? What sorts of inferential questions can we ask and answer?

We do still need to make *some* assumptions (of course), but they can be weaker than what we're used to. For example, instead of normality, which is a strong assumption, we might assume that the true data distribution is merely symmetric.

For comparing two samples, rather than assuming that both come from normally-distributed populations with possibly different means, we might assume that their distributions are the same (without specifying what it is) but with a **shift** in location:

```r
library(tidyverse)

N <- 1e+6

components <- sample(1:3,size = N,replace = TRUE,
                     prob = c(0.3,0.5,0.2))
mus <- c(0,10,3)
sds <- sqrt(c(0.2,1,3))

samples <- rnorm(
  n = N,
  mean = mus[components],
  sd = sds[components]
)

tibble(
  weird_1 = samples,
  weird_2 = samples + 2
) %>%
  gather(key = "dist", value = "value") %>%
  ggpubr::ggdensity(x = "value", color = "dist", fill = "dist",
                    palette = "npg", alpha = 0.25)
```

{{< figure src="1-arbitrary_dist.png" title="Two arbitrary distributions with the same shape but different means." numbered="true" lightbox="true" >}}



[^1]: methods will be good in a wider range of applications.