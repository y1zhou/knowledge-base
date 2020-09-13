---
title: "ARMA Model"
date: 2020-09-12T20:31:18-04:00
summary: "The mean, variance, ACF and PACF of " # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "time-series-arma-model"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 23 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

## Backshift operator

The `backshift operator` $B$ shifts the time index to the previous one[^backshift-operator]:

$$
BX_t = X_{t-1}, B^2 X_t = X_{t-2}, \cdots
$$

In general, $B^k X_t = X_{t-k}$ represents shifting $X_t$ $k$ units back in time. The backshift operator doesn't operate on coefficients because they are fixed quantities that do not move in time.

[^backshift-operator]: Note that $B$ is an operator, so $X_t B$ is invalid.

### AR models

For an AR(1) model with $\mu = 0$,

$$
\begin{gathered}
    X_t = \phi X_{t-1} + Z_t \\\\
    X_t - \phi X_{t-1} = Z_t \Rightarrow X_t - \phi BX_t = Z_t \\\\
    (1 - \phi B)X_t = Z_t
\end{gathered}
$$

For an AR(p) model with $\mu = 0$:

$$
\begin{gathered}
    X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \cdots + \phi_p X_{t-p} + Z_t \\\\
    X_t = \phi_1 BX_t + \phi_2 B^2 X_t + \cdots + \phi_p B^p X_t + Z_t \\\\
    X_t - \phi_1 BX_t - \phi_2 B^2 X_t - \cdots - \phi_p B^p X_t = Z_t \\\\
    (1 - \phi_1 B - \phi_2 B^2 - \cdots - \phi_p B^p)X_t = Z_t
\end{gathered}
$$

The general form for an `AR polynomial` is:

$$
\phi(B) = 1 - \phi_1 B - \cdots - \phi_p B^p
$$

So writing an AR(p) model using the AR polynomial:

$$
\phi(B) X_t = \delta + Z_t
$$

### MA models

For an MA(q) model with $\mu = 0$,

$$
\begin{gathered}
    X_t = Z_t + \theta_1 Z_{t-1} + \theta_2 Z_{t-2} + \cdots + \theta_q Z_{t-q} \\\\
    X_t = Z_t + \theta_1 BZ_t + \theta_2 B^2 Z_t + \cdots + \theta_q B^q Z_t \\\\
    X_t = (1 + \theta_1 B + \theta_2 B^2 + \cdots + \theta_q B^q)Z_t
\end{gathered}
$$

Using the `MA polynomial`, this can be expressed as

$$
X_t = \theta(B)Z_t \text{ where } \theta(B) = 1 + \theta_1 B + \cdots + \theta_q B^q
$$

The condition for invertibiliity of an MA(q) series is the following equation

$$
1 - \theta_1 z - \theta_2 z^2 - \cdots - \theta_q z^q = 0
$$

has solutions for $z$ that fall outside the unit circle.

## ARMA

As we've seen in some of the ACF and PACF plots, it's often hard to distinguish between "cut off" and "tail off" patterns. In an `ARMA` model, we consider both functions to tail off. ARMA is more flexible as it combines MA and AR models into one, but the downside is more parameters are added to the model. In general, large number of parameters reduces the efficiency in estimation.

### ARMA(1, 1)

In an `ARMA(1, 1)` model, the LHS looks like an AR(1) model and the RHS resembles an MA(1) model:

$$
X_t - \phi X_{t-1} = \delta + Z_t + \theta Z_{t-1}
$$

where $\\{Z_t\\} \sim WN(0, \sigma^2)$. We can express it using AR and MA polynomials:

$$
(1 - \phi B)X_t = \delta + (1 + \theta B)Z_t
$$

We can also drop the intercept term and incorporate the mean $\mu$:

$$
(1 - \phi B)(X_t - \mu) = (1 + \theta B)Z_t
$$
