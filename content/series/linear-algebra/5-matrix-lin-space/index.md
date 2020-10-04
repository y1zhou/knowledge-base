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
    \mathcal{C}(\boldsymbol{A}) &= \mathcal{L}(\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n) \\\\
    &= \left\\{ \sum\_{j=1}^n c_j \boldsymbol{a}_j, \quad c_j \in \mathbb{R} \right\\} \\\\
    &= \\{ \boldsymbol{Ac}, \quad \boldsymbol{c} \in \mathbb{R}^n \\}, \quad \boldsymbol{c} = (c_1, \cdots, c_n)^\prime
\end{aligned}
$$

This means {{<hl>}}for any $\boldsymbol{y} \in \mathcal{C}(\boldsymbol{A})$, there exists $\boldsymbol{x} \in \mathbb{R}^n$ such that $\boldsymbol{y} = \boldsymbol{Ax}$.{{</hl>}}

Similarly, the `row space` of $\boldsymbol{A}$ is

$$
\begin{aligned}
    \mathcal{R}(\boldsymbol{A}) &= \mathcal{L}(\boldsymbol{\alpha}_1^\prime, \cdots, \boldsymbol{\alpha}_m^\prime) \\\\
    &= \left\\{ \sum\_{i=1}^m c_i \boldsymbol{\alpha}_j^\prime, \quad c_i \in \mathbb{R} \right\\} \\\\
    &= \\{ \boldsymbol{c}^\prime \boldsymbol{A}, \quad \boldsymbol{c} \in \mathbb{R}^m \\}
\end{aligned}
$$

Note that

$$
\begin{gathered}
    \boldsymbol{y} \in \mathcal{C}(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{y}^\prime \in \mathcal{R}(\boldsymbol{A}^\prime) \\\\
    \boldsymbol{x}^\prime \in \mathcal{R}(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{x} \in \mathcal{C}(\boldsymbol{A}^\prime)
\end{gathered}
$$

### Properties

1. $dim(\mathcal{C}(\boldsymbol{A})) \leq \min(n, m)$. The dimension of the column space of $\boldsymbol{A}$ is called the `rank` of $\boldsymbol{A}$, which is written as $rank(\boldsymbol{A}) \triangleq r(\boldsymbol{A})$.
2. $dim(\mathcal{R}(\boldsymbol{A})) = dim(\mathcal{C}(\boldsymbol{A})) = r(\boldsymbol{A})$, which implies the number of LIN columns is equal to the number of LIN rows.

    **Proof**: Show that the row rank is equal to the column rank. Let $r$ denote the row rank, that is we have $r$ linearly independent rows. This means $\mathcal{R}(\boldsymbol{A})$ has a basis $\\{\boldsymbol{x}_1^\prime, \cdots, \boldsymbol{x}_r^\prime\\}$. We will show that $\\{\boldsymbol{Ax}_1, \cdots, \boldsymbol{Ax}_r\\}$ is linearly independent.

    We show this because each vector $\boldsymbol{Ax}_i \in \mathcal{C}(\boldsymbol{A})$. If we can show the set is linearly independent, then the dimension of the column space has to be at least $r$. To show this, we set

    $$
    \sum_{i=1}^r c_i \boldsymbol{Ax}_i = \boldsymbol{0}
    $$

    where the $c_i$'s are scalars. If we take $\boldsymbol{A}$ outside,

    $$
    \boldsymbol{A} \left( \sum_{i=1}^r c_i \boldsymbol{x}_i \right) = \boldsymbol{0}
    $$

    Let $\boldsymbol{v} = \sum_{i=1}^r c_i \boldsymbol{x}_i$. If $\boldsymbol{Av} = \boldsymbol{0}$, then $\boldsymbol{v}^\prime$ has to be orthogonal to each row of $\boldsymbol{A}$, which means $\boldsymbol{v}^\prime \perp \mathcal{R}(\boldsymbol{A})$. Note that $\boldsymbol{v}^\prime$ is a linear combination of the rows of $\boldsymbol{A}$, so $\boldsymbol{v}^\prime \in \mathcal{R}(\boldsymbol{A})$. The only vector that belongs to a subspace and is orthogonal to that subspace is the zero vector.

    If $\boldsymbol{v} = \boldsymbol{0}$, then $c_i = 0$ for all $i$ because the $\boldsymbol{x}_i$'s are linearly independent. Thus, the $\boldsymbol{Ax}_i$'s are linearly independent. The column space of $\boldsymbol{A}$ has a set of $r$ linearly independent vectors, so $dim(\mathcal{C}(\boldsymbol{A})) \geq r$.

    We can repeat the above argument with $c \triangleq dim(\mathcal{C}(\boldsymbol{A}))$, and show that $dim(\mathcal{R}(\boldsymbol{A})) \geq c$. By showing both sides of the inequality, we arrive at the conclusion that $r=c$, i.e. the dimension of the column space is the same as the dimension of the row space, and that's what we call the rank of $\boldsymbol{A}$.

3. If $r(\boldsymbol{A}) = m$, then we say $\boldsymbol{A}$ has full row rank[^full-row-rank]. If $r(\boldsymbol{A}) = n$, then we say $\boldsymbol{A}$ has full column rank. Either case, we say $\boldsymbol{A}$ has `full rank`. If a matrix is not of full rank, it is `rank-deficient`.

    [^full-row-rank]: It's more likely a short, wide matrix. Similarly, full column rank matrices are likely to be narrow and tall.

4. If $r(\boldsymbol{A}) = r$, then there exists matrix $\boldsymbol{B}: m \times r$, $r(\boldsymbol{B}) = r$ (full column rank) and matrix $\boldsymbol{T}: r \times n$, $r(\boldsymbol{T}) = r$ (full row rank), such that $\boldsymbol{A} = \boldsymbol{BT}$.

    Any matrix can be written as a product of two full rank matrices. $\mathcal{C}(\boldsymbol{A})$ has $r$ basis vectors, and we can put them into an $m \times r$ matrix $\boldsymbol{B}$. Each column of $\boldsymbol{A}$, $\boldsymbol{a}_i$, is in $\mathcal{C}(\boldsymbol{A})$, so $\boldsymbol{a}_i$ can be expressed as $\boldsymbol{a}_i = \boldsymbol{By}_i$ where the $\boldsymbol{y}_i$'s are $r \times 1$ vectors, and we can put them together into an $r \times n$ matrix $\boldsymbol{T}$.

5. $r(\boldsymbol{AB}) \leq \min(r(\boldsymbol{A}), r(\boldsymbol{B}))$ for any $\boldsymbol{A}$, $\boldsymbol{B}$ matrices.

    **Proof**: We first show that $r(\boldsymbol{AB}) \leq r(\boldsymbol{A})$ from the fact that $\mathcal{C}(\boldsymbol{AB}) \subset \mathcal{C}(\boldsymbol{A})$. Note that

    $$
    \boldsymbol{AB} = [\boldsymbol{Ab}_1, \cdots, \boldsymbol{Ab}_q]
    $$

    where $\boldsymbol{B}$ has $q$ columns. Each column of $\boldsymbol{AB}$ has a form of $\boldsymbol{Ab}_i$, which is in $\mathcal{C}(\boldsymbol{A})$. In other words, any vector that lives in the column space of $\boldsymbol{AB}$ is in the column space of $\boldsymbol{A}$.

    Similarly, we can show $r(\boldsymbol{AB}) \leq r(\boldsymbol{B})$ from the fact that $\mathcal{R}(\boldsymbol{AB}) \subset \mathcal{R}(\boldsymbol{B})$[^column-space-of-ab].

    [^column-space-of-ab]: $\mathcal{C}(\boldsymbol{AB})$ and $\mathcal{C}(\boldsymbol{B})$ can't even be compared usless $m = n$. However their dimensions could be discussed.

6. If an $m \times n$ matrix $\boldsymbol{A}$ with rank $r$ is written as $\boldsymbol{A} = \boldsymbol{BT}$ with $\boldsymbol{B}: m \times r$ and $\boldsymbol{T}: r \times n$, then $r(\boldsymbol{B}) = r(\boldsymbol{T}) = r$.