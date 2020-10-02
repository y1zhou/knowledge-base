---
title: "Linear Space of Matrices"
date: 2020-09-30T13:23:18-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "linear-algebra-matrix-linear-space"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: true

weight: 50 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

In the previous two sections (vector space and geometrical considerations) we've been focusing on vector operations and properties. Now we're ready to come back to matrices. Associated with any matrix is a very important characteristic called the rank. There are several (consistent) ways of defining the rank. The most fundamental of these is in terms of the dimension of a linear space.

## Definitions

Recall a matrix $\boldsymbol{A}: m \times n$ has $\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n \in \mathbb{R}^m$ columns and $\boldsymbol{\alpha}_1^\prime, \cdots, \boldsymbol{\alpha}_m^\prime \in \left( \mathbb{R}^n \right)^\prime$ rows. The `column space` of $\boldsymbol{A}$ contains all linear combinations of the columns:

$$
\begin{aligned}
    C(\boldsymbol{A}) &= \mathcal{L}(\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n) \\\\
    &= \left\\{ \sum\_{j=1}^n c_j \boldsymbol{a}_j, \quad c_j \in \mathbb{R} \right\\} \\\\
    &= \\{ \boldsymbol{Ac}, \quad \boldsymbol{c} \in \mathbb{R}^n \\}, \quad \boldsymbol{c} = (c_1, \cdots, c_n)^\prime
\end{aligned}
$$

This means {{<hl>}}for any $\boldsymbol{y} \in C(\boldsymbol{A})$, there exists $\boldsymbol{x} \in \mathbb{R}^n$ such that $\boldsymbol{y} = \boldsymbol{Ax}$.{{</hl>}}

Similarly, the `row space` of $\boldsymbol{A}$ is

$$
\begin{aligned}
    R(\boldsymbol{A}) &= \mathcal{L}(\boldsymbol{\alpha}_1^\prime, \cdots, \boldsymbol{\alpha}_m^\prime) \\\\
    &= \left\\{ \sum\_{i=1}^m c_i \boldsymbol{\alpha}_j^\prime, \quad c_i \in \mathbb{R} \right\\} \\\\
    &= \\{ \boldsymbol{c}^\prime \boldsymbol{A}, \quad \boldsymbol{c} \in \mathbb{R}^m \\}
\end{aligned}
$$

Note that

$$
\begin{gathered}
    \boldsymbol{y} \in C(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{y}^\prime \in R(\boldsymbol{A}^\prime) \\\\
    \boldsymbol{x}^\prime \in R(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{x} \in C(\boldsymbol{A}^\prime)
\end{gathered}
$$
