---
title: "Matrix Trace"
date: 2020-10-07T13:07:16-04:00
summary: "Such a simple concept with so many properties and applications!" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix", "Operation", "Statistics"] # list of related tags

slug: "linear-algebra-matrix-trace"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 60 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Night sky with stars." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "SGX0EAkJyqQ" # Unsplash ID of the picture
---

The `trace` of a square matrix $\boldsymbol{A}: n \times n$ is the sum of the diagonals:

$$
tr(\boldsymbol{A}) = \sum_{i=1}^n a_{ii}
$$

The trace is related to many important features of matrices, and we're going to discuss some basic properties first.

## Properties

-   **Trace is commutative**, i.e. $tr(\boldsymbol{AB}) = tr(\boldsymbol{BA})$. Note that if $\boldsymbol{A}: m \times n$ and $\boldsymbol{B}: n \times m$, $\boldsymbol{AB}$ would be $m \times m$ and $\boldsymbol{BA}$ would be $n \times n$.

$$
\begin{aligned}
    LHS &= \sum_{i=1}^m (\boldsymbol{AB})_{ii} \\\\
    &= \sum_{i=1}^m \left( \sum_{i=1}^n a_{ij}b_{ji} \right) \\\\
    &= \sum_{j=1}^n \sum_{i=1}^m b_{ji}a_{ij} \\\\
    &= \sum_{j=1}^n (\boldsymbol{BA})_{jj} = tr(\boldsymbol{BA})
\end{aligned}
$$

For more than two matrices, note that $tr(\boldsymbol{ABC}) \neq tr(\boldsymbol{CBA})$. Instead,

$$
tr(\boldsymbol{ABC}) = tr(\boldsymbol{A}(\boldsymbol{BC})) = tr(\boldsymbol{BCA})
$$

Similarly, we have

$$
tr(\boldsymbol{ABCDE}) = tr(\boldsymbol{CDEAB}) = tr(\boldsymbol{EABCD}) = \cdots
$$

-   $tr(\boldsymbol{A}^\prime \boldsymbol{A}) = tr(\boldsymbol{AA}^\prime) \geq 0$.

This is equal to the sum of all entries squared.

-   $tr(\boldsymbol{A}^\prime\boldsymbol{A}) = 0$ is equivalent to $\boldsymbol{A} = 0$.

Showing $\Leftarrow$ is trivial. For $\Rightarrow$, the trace being zero means that the sum of all entries squared of $\boldsymbol{A}$ is zero, thus $a_{ij} = 0$ for all $i, j$.

We can also use this fact to show that $tr(\boldsymbol{A}^\prime \boldsymbol{A}) = 0 \Leftrightarrow \boldsymbol{A}^\prime\boldsymbol{A} = 0$.

-   $\boldsymbol{A}^\prime \boldsymbol{AB} = \boldsymbol{A}^\prime \boldsymbol{AC} \Leftrightarrow \boldsymbol{AB} = \boldsymbol{AC}$.

This one is related to future statistical analyses (e.g. linear regression). Again $\Leftarrow$ is obvious. To prove $\Rightarrow$, we will show $\boldsymbol{AB} = \boldsymbol{AC}$ by seeing that $tr((\boldsymbol{AB} - \boldsymbol{AC})^\prime (\boldsymbol{AB} - \boldsymbol{AC})) = 0$.

We start with $\boldsymbol{A}^\prime \boldsymbol{AB} - \boldsymbol{A}^\prime \boldsymbol{AC} = 0$.

$$
\begin{gathered}
    \boldsymbol{A}^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    (\boldsymbol{B}^\prime - \boldsymbol{C}^\prime) \boldsymbol{A}^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    (\boldsymbol{AB} - \boldsymbol{AC})^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    \boldsymbol{AB} - \boldsymbol{AC} = 0
\end{gathered}
$$

We will use this fact to show that $\mathcal{C}(\boldsymbol{A}) = \mathcal{C}(\boldsymbol{AA}^\prime)$.

### Recap on column spaces

Say $\boldsymbol{A}$ is an $m \times n$ matrix. If vector $\boldsymbol{y} \in \mathcal{C}(\boldsymbol{A})$, it means $\boldsymbol{y} = \boldsymbol{Ax}$ for some $\boldsymbol{x}$. Similarly if $\boldsymbol{y} \in \mathcal{C}(\boldsymbol{A})^\perp$, then $\boldsymbol{y} \perp \mathcal{C}(\boldsymbol{A})$. Since $\boldsymbol{y}$ is orthogonal to each column of $\boldsymbol{A}$, we have

$$
\boldsymbol{y}^\prime \boldsymbol{a}_j = 0 \quad \forall j = 1, \cdots, n \Rightarrow \boldsymbol{y}^\prime \boldsymbol{A} = \boldsymbol{0}^\prime
$$

Similarly, if $\boldsymbol{y}^\prime \in \mathcal{R}(\boldsymbol{A})$, then $\boldsymbol{y}^\prime = \boldsymbol{x}^\prime \boldsymbol{A}$ for some $\boldsymbol{x}$, and

$$
\boldsymbol{y}^\prime \in \underbrace{\mathcal{R}(\boldsymbol{A})^\perp}\_{=N(\boldsymbol{A})} \Leftrightarrow \boldsymbol{Ay} = \boldsymbol{0}_{m \times 1}
$$

-   Suppose $\mathcal{W}, \mathcal{V}$ are subspaces in $\mathbb{R}^n$, and $\mathcal{W} \subset \mathcal{V}$. Then, $\mathcal{V}^\perp \subset \mathcal{W}^\perp$.

We can express the subspaces as spans:

$$
\begin{gathered}
    \mathcal{W} = \mathcal{L}(\boldsymbol{u}_1, \cdots, \boldsymbol{u}_k) \\\\
    \mathcal{V} = \mathcal{L}(\boldsymbol{u}_1, \cdots, \boldsymbol{u}_k, \boldsymbol{u}\_{k+1}, \cdots, \boldsymbol{u}_m)
\end{gathered}
$$

For any $\boldsymbol{y} \in \mathcal{V}^\perp$, we have

$$
\boldsymbol{y} \perp \mathcal{V} \Leftrightarrow \boldsymbol{y} \perp \boldsymbol{u}_i, i = 1, \cdots, m
$$

Which means $\boldsymbol{y} \perp \boldsymbol{u}_i, i = 1, \cdots, k$, thus $\boldsymbol{y} \perp \mathcal{W}$ and $\boldsymbol{y} \in \mathcal{W}^\perp$.

### Proof

We're now ready to prove $\mathcal{C}(\boldsymbol{A}^\prime) = \mathcal{C}(\boldsymbol{A}^\prime \boldsymbol{A})$. Usually is not that easy to prove set equality directly. Often we prove the LHS is a subset of the RHS and vice versa.

We first show that $\mathcal{C}(\boldsymbol{A}^\prime \boldsymbol{A}) \subset \mathcal{C}(\boldsymbol{A}^\prime)$. This is obvious because $\mathcal{C}(\boldsymbol{AB}) \subset \mathcal{C}(\boldsymbol{A})$.

To show that $\mathcal{C}(\boldsymbol{A}^\prime) \subset \mathcal{C}(\boldsymbol{A}^\prime \boldsymbol{A})$, we will show $\mathcal{C}(\boldsymbol{A}^\prime \boldsymbol{A})^\perp \subset \mathcal{C}(\boldsymbol{A}^\prime)^\perp$ instead. For all $\boldsymbol{y} \in \mathcal{C}(\boldsymbol{A}^\prime \boldsymbol{A})^\perp$,

$$
\begin{gathered}
    \quad \boldsymbol{y}^\prime \boldsymbol{A}^\prime \boldsymbol{A} = \boldsymbol{0}^\prime
    \xrightarrow{transpose} \boldsymbol{A}^\prime \boldsymbol{Ay} = \boldsymbol{0}\_{n \times 1} \\\\
    \Rightarrow \boldsymbol{Ay} = \boldsymbol{0}\_{m \times 1} \xrightarrow{transpose} \boldsymbol{y}^\prime \boldsymbol{A} = \boldsymbol{0}^\prime \\\\
    \boldsymbol{y} \in \mathcal{C}(\boldsymbol{A}^\prime)^\perp
\end{gathered}
$$

Similarly, we can show that

$$
\begin{gathered}
    \mathcal{C}(\boldsymbol{A}) = \mathcal{C}(\boldsymbol{AA}^\prime) \\\\
    \mathcal{R}(\boldsymbol{A}) = \mathcal{R}(\boldsymbol{A}^\prime \boldsymbol{A}) \\\\
    r(\boldsymbol{A}) = r(\boldsymbol{AA}^\prime) = r(\boldsymbol{A}^\prime\boldsymbol{A})
\end{gathered}
$$

## One-way ANOVA example

Suppose matrix $\boldsymbol{A}$ is given by

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 1 & 0 \\\\
    1 & 1 & 0 \\\\
    1 & 0 & 1 \\\\
    1 & 0 & 1
\end{pmatrix}
$$

We can easily see that the second the third columns are linearly independent, but the first column is the sum of the other two columns, so $r(\boldsymbol{A}) = 2$. The column space of $\boldsymbol{A}$ is a subspace of $\mathbb{R}^4$. Suppose the data we collected is

$$
\boldsymbol{y} = \begin{pmatrix}
    2 \\\\ 1 \\\\ 3 \\\\ 5
\end{pmatrix}
$$

In ANOVA, we're going to need two more matrices:

$$
\boldsymbol{A}^\prime \boldsymbol{A} = \begin{pmatrix}
    4 & 2 & 2 \\\\
    2 & 2 & 0 \\\\
    2 & 0 & 2
\end{pmatrix}, \quad \boldsymbol{A}^\prime \boldsymbol{y} = \begin{pmatrix}
    11 \\\\ 3 \\\\ 8
\end{pmatrix}
$$

and we need to solve for $\boldsymbol{\beta} \in \mathbb{R}^3$:

$$
\boldsymbol{A}^\prime \boldsymbol{A\beta} = \boldsymbol{A}^\prime \boldsymbol{y}
$$

We're interested in how to solve this, if the solution is unique, or does a solution even exist. The next chapter is dedicated to linear systems.
