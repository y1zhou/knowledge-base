---
title: "Matrix Inverse"
date: 2020-10-12T12:42:03-04:00
summary: "... for a nonsingular matrix." # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra"] # list of related tags

slug: "linear-algebra-matrix-inverse"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 70 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

## Linear system

We have matrix $\boldsymbol{A}: m \times n$, and vectors $\boldsymbol{b}: m \times 1$ and $\boldsymbol{x}: n \times 1$ where $\boldsymbol{b}$ is given and $\boldsymbol{x}$ is unknown. We use the equation

$$
\boldsymbol{Ax} = \boldsymbol{b}
$$

We're interested in when does a solution $\boldsymbol{x}$ exist, and how do we find it.

A linear system is `consistent` when there exists at least one $\boldsymbol{x}$ satisfying the equation. We can see that the following statements are equivalent:

1. A linear system is consistent.
2. $\boldsymbol{b} \in \mathcal{C}(\boldsymbol{A})$.
3. $\mathcal{C}([\boldsymbol{A}, \boldsymbol{b}]) = \mathcal{C}(\boldsymbol{A})$.
4. $r([\boldsymbol{A}, \boldsymbol{b}]) = r(\boldsymbol{A})$.

If $\boldsymbol{A}$ has full row rank $m$, then $\mathcal{C}(\boldsymbol{A}) = \mathbb{R}^m$, and $\boldsymbol{Ax} = \boldsymbol{b}$ is consistent for any $\boldsymbol{b}$.

The next fact is important in statistical contexts. The linear system $\boldsymbol{A}'\boldsymbol{Ax} = \boldsymbol{A}'\boldsymbol{b}$ is consistent for any $\boldsymbol{A}$, $\boldsymbol{b}$. This is because $\boldsymbol{A}'\boldsymbol{b} \in \mathcal{C}(\boldsymbol{A}') = \mathcal{C}(\boldsymbol{A}'\boldsymbol{A})$.

If $\boldsymbol{A}$ has full column rank, $\boldsymbol{A}'\boldsymbol{A}$ is nonsingular, and $\boldsymbol{x}$ is unique. This is an assumption in our regression analyses, which is there's no redundancy in our predictor variables.

## Matrix inverse

The inverse of a matrix only makes sense when the matrix is square. Before talking about it, we introduce two other concepts for non-square matrices.

### Left and right inverses

If $\boldsymbol{A}: m \times n$, the `right inverse` $\boldsymbol{R}$ is the $n \times m$ matrix such that

$$
\boldsymbol{AR} = \boldsymbol{I}_m
$$

The right inverse exists when $\boldsymbol{A}$ has full row rank, i.e. $r(\boldsymbol{A}) = m$ ($m \leq n$). This is because if $r(\boldsymbol{A}) = m$, then $\boldsymbol{Ax} = \boldsymbol{b}$ is consistent for any $\boldsymbol{b}$. Let

$$
\boldsymbol{b}_1 = \begin{pmatrix}
    1 \\\\ 0 \\\\ 0 \\\\ \vdots \\\\ 0
\end{pmatrix}\_{m \times 1},
\boldsymbol{b}_2 = \begin{pmatrix}
    0 \\\\ 1 \\\\ 0 \\\\ \vdots \\\\ 0
\end{pmatrix}, \cdots,
\boldsymbol{b}_m = \begin{pmatrix}
    0 \\\\ 0 \\\\ 0 \\\\ \vdots \\\\ 1
\end{pmatrix}
$$

The solution for $\boldsymbol{Ab}_1 = \boldsymbol{x}$ is $\boldsymbol{x}_1$, and similarly we have $\boldsymbol{x}_2, \cdots, \boldsymbol{x}_m$. If we put $\boldsymbol{x}_1, \cdots, \boldsymbol{x}_m$ into an $n \times m$ matrix

$$
\boldsymbol{R} = [\boldsymbol{x}_1, \cdots, \boldsymbol{x}_m],
$$

we have $\boldsymbol{AR} = \boldsymbol{A}[\boldsymbol{x}_1, \cdots, \boldsymbol{x}_m] = \boldsymbol{I}_m$. From the proof we can also see that $\boldsymbol{R}$ is not unique.

Similarly, the `left inverse` $\boldsymbol{L}: n \times m$ is the matrix that satisfies $\boldsymbol{LA} = \boldsymbol{I}_n$. $\boldsymbol{L}$ exists when $\boldsymbol{A}$ has full column rank, i.e. $r(\boldsymbol{A}) = n (m \geq n)$.

#### Inverse

If $\boldsymbol{A}$ is non-singular, then its right inverse and left inverse exist, and they are the same, which we call the `inverse` of $\boldsymbol{A}$ denoted by $\boldsymbol{A}^{-1}$. $\boldsymbol{A}^{-1}$ is unique. To prove this, note that

$$
\boldsymbol{L} = \boldsymbol{LI}_n = \boldsymbol{LAR} = \boldsymbol{I}_n \boldsymbol{R} = \boldsymbol{R}
$$

#### Properties

If $\boldsymbol{A}$ has full column rank $n$, and there exists $\boldsymbol{L}$ such that $\boldsymbol{LA} = \boldsymbol{I}_n$ ($n \leq m$), then

1. $\boldsymbol{B} = \boldsymbol{C} \Longleftrightarrow \boldsymbol{AB} = \boldsymbol{AC}$

    This is because

    $$
    \boldsymbol{AB} = \boldsymbol{AC} \Rightarrow \boldsymbol{LAB} = \boldsymbol{LAC} \Rightarrow \boldsymbol{I}_n\boldsymbol{B} = \boldsymbol{I}_n \boldsymbol{C}
    $$

2. For any $\boldsymbol{B}$, $r(\boldsymbol{AB}) = r(\boldsymbol{B})$. Pre-multiplying a full column rank matrix does not change the rank.
