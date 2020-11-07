---
title: "Quadratic Form"
date: 2020-11-04T16:45:31-05:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra"] # list of related tags

slug: "linear-algebra-quadratic-form"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 110 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Ship breaking down icebergs." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "5zFxCK_p6So" # Unsplash ID of the picture
---

The quadratic form of a matrix is a function of some vector. We're going to discuss several cases, starting with the simplest case of a linear form.

## Definitions

### Linear form

The linear form of $\boldsymbol{x} = (x_1, \cdots, x_n)'$ is

$$
\begin{gathered}
    \boldsymbol{a}'\boldsymbol{x} = a_1 x_1 + a_2 x_2 + \cdots + a_n x_n \\\\
    \boldsymbol{Ax} = \begin{pmatrix}
        \boldsymbol{\alpha}_1'\boldsymbol{x} \\\\
        \vdots \\\\
        \boldsymbol{\alpha}_m'\boldsymbol{x}
    \end{pmatrix}
\end{gathered}
$$

where $\boldsymbol{\alpha}_i'$ denotes the $i$-th row of $\boldsymbol{A}$. Each row of $\boldsymbol{Ax}$ is a linear combination of the elements of $\boldsymbol{x}$, which is why it is called a `linear functional` or `linear form`. Note that $\boldsymbol{A}$ is a given matrix, and $\boldsymbol{x}$ is just a vector with known dimensions[^understanding-x].

[^understanding-x]: Just like when we first learned about solving $x-1=0$, where $x$ denotes some unknown number.

### Bilinear form

Taking this a step further, we now have two vectors $\boldsymbol{x}$ and $\boldsymbol{y}$:

$$
\boldsymbol{x}'\boldsymbol{Ay}, \quad \boldsymbol{A}: m \times n, \boldsymbol{x}: m \times 1, \boldsymbol{y}: n \times 1
$$

We're multiplying a row vector to a matrix, then two a column vector, so the outcome is going to be a scalar value. Again, only $\boldsymbol{A}$ is a known matrix, and the vectors $\boldsymbol{x}$ and $\boldsymbol{y}$ are just unknowns.

$$
\begin{aligned}
    \boldsymbol{x}'\boldsymbol{Ay} &= (x_1, \cdots, x_m) \begin{pmatrix}
        a_{11} & \cdots & a_{1n} \\\\
        \vdots & \ddots & \vdots \\\\
        a_{m1} & \cdots & a_{mn}
    \end{pmatrix}
    \begin{pmatrix}
        y_1 \\\\ \vdots \\\\ y_n
    \end{pmatrix} \\\\
    &= \sum_{j=1}^n \sum_{i=1}^m a_{ij}x_i y_j
\end{aligned}
$$

For example, let

$$
\boldsymbol{A} = \begin{pmatrix}
    2 & 1 & 3 \\\\
    -1 & 0 & 1
\end{pmatrix}, \quad \boldsymbol{x}'\boldsymbol{Ay} = 2x_1y_1 + x_1y_2 + 3x_1y_3 - x_2y_1 + x_2y_3
$$

For some fixed $\boldsymbol{x}$, it becomes a linear function of $\boldsymbol{y}$. Similarly if $\boldsymbol{y}$ is fixed, then it becomes a linear function of $\boldsymbol{x}$.

### Quadratic form

The difference the bilinear form and the quadratic form is that now $\boldsymbol{x}$ and $\boldsymbol{y}$ are equal to each other, and thus $\boldsymbol{A}$ must be a square matrix.

Any $n \times n$ real matrix $\boldsymbol{A}$ determines a quadratic form in $n$ variables $\boldsymbol{x} = (x_1, \cdots, x_n)'$ by the formula

$$
\boldsymbol{x}'\boldsymbol{Ax} = \sum_{j=1}^n \sum_{i=1}^n a_{ij}x_i x_j = \sum_{i=1}^n a_{ii}x_i^2 + \sum_{i \neq j} a_{ij}x_i x_j
$$

For example, let

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 2 \\\\
    3 & 4
\end{pmatrix}
$$

The quadratic form is then

$$
\boldsymbol{x}'\boldsymbol{Ax} = x_1^2 + 4x_2^2 + 2x_1x_2 + 3x_2x_1 = x_1^2 + 4x_2^2 + 5x_1x_2
$$

So why are we doing this? The point is that investigating properties of this function reveals properties about the matrix itself. Knowing the function $x_1^2 + 4x_2^2 + 5x_1x_2$ in the space $\mathbb{R}^2$ tells us a lot about $\boldsymbol{A}$.

### Notes on the quadratic form

1. $\boldsymbol{x}'\boldsymbol{Ax} = \boldsymbol{x}'\boldsymbol{A}'\boldsymbol{x}$. The LHS $\boldsymbol{x}'\boldsymbol{Ax}$ is a scalar and thus symmetric, so

    $$
    \left(\boldsymbol{x}'\boldsymbol{Ax}\right)' = \boldsymbol{x}'\boldsymbol{A}'\boldsymbol{x} = RHS
    $$

    This implies that as far as quadratic forms are concerned, we can confine to symmetric matrices, i.e. for any $n \times n$ matrix $\boldsymbol{A}$, there exists a symmetric matrix $\boldsymbol{B}$ that has exactly the same quadratic forms.

    $$
    \boldsymbol{B} = \frac{\boldsymbol{A} + \boldsymbol{A}'}{2}, \quad \boldsymbol{x}'\boldsymbol{Bx} = \frac{1}{2}\boldsymbol{x}'\boldsymbol{Ax} + \frac{1}{2} \boldsymbol{x}'\boldsymbol{A}'\boldsymbol{x} = \boldsymbol{x}'\boldsymbol{Ax}
    $$

2. If $\boldsymbol{A}$ and $\boldsymbol{B}$ are symmetric matrices, then

    $$
    \boldsymbol{x}'\boldsymbol{Ax} = \boldsymbol{x}'\boldsymbol{Bx} \Longleftrightarrow \boldsymbol{A} = \boldsymbol{B}
    $$

## Categorizing quadratic forms

Now we're ready to talk about categorizing matrices, or equivalently categorizing quadratic forms. The concepts we're going to introduce come from a desire for notions about the signs of matrices. In elementary maths, we have $3 > 0$ and $-3 < 0$, but for two matrices:

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 2 \\\\
    2 & 3
\end{pmatrix} \quad vs. \quad
\boldsymbol{B} = \begin{pmatrix}
    2 & 1 \\\\
    1 & 2
\end{pmatrix},
$$

we can't easily say one is greater than the other. The notion of definiteness will help us with this comparison. Conceptually, we have

$$
\boldsymbol{A} - \boldsymbol{B} = \begin{pmatrix}
    -1 & 1 \\\\
    1 & 1
\end{pmatrix},
$$

and then we determine the definiteness of the matrix $\boldsymbol{A} - \boldsymbol{B}$ as a whole, and say "$\boldsymbol{A}$ is larger than $\boldsymbol{B}$" or the other way around.

### Non-negative definite

The definition of `non-negative definite` (NND) matrices applies to quadratic forms as well as to matrices directly. A quadratic form $\boldsymbol{x}'\boldsymbol{Ax}$ is NND (equivalently $\boldsymbol{A}$ is NND) if $\boldsymbol{x}'\boldsymbol{Ax} \geq 0$ for any $\boldsymbol{x} \in \mathbb{R}^n$.

Clearly, $\boldsymbol{x} = 0$ yields $\boldsymbol{x}'\boldsymbol{Ax} = 0$. Two distinct cases for NND are:

1. $\boldsymbol{x} = \boldsymbol{0}$ is the **only** vector that yields zero. This case is called `positive definite` (PD).
2. There are some $\boldsymbol{x} \neq \boldsymbol{0}$ that yield $\boldsymbol{x}'\boldsymbol{Ax} = 0$. This case is called `positive semi-definite` (PSD).

#### Examples

-   The identity matrix is positive definite.

    $$
    \boldsymbol{x}'\boldsymbol{Ix} = x_1^2 + x_2^2 + \cdots + x_n^2 \geq 0
    $$

-   The matrix of ones is positive semi-definite.

    $$
    \begin{aligned}
        \boldsymbol{x}'\boldsymbol{Jx} &= x_1^2 + \cdots + x_n^2 + 2(x_1 x_2 + \cdots + x_{n-1}x_n) \\\\
        &= (x_1 + \cdots + x_n)^2 \geq 0 \quad \Rightarrow NND
    \end{aligned}
    $$

    This can be zero if $x_1 + \cdots + x_n = 0$, so $\boldsymbol{J}$ is not positive-definite, and is thus positive semi-definite.

### Non-positive definite

On the opposite side of NND matrices, we have non-positive definite (NPD) matrices. We say a quadratic form $\boldsymbol{x}'\boldsymbol{Ax}$ is `non-positive definite` if $\boldsymbol{x}'\boldsymbol{Ax} \leq 0$. The two possible cases are:

1. `Negative definite`: $\boldsymbol{x} = \boldsymbol{0}$ is the only vector to make $\boldsymbol{x}'\boldsymbol{Ax} = 0$.
2. `Negative semi-definite`: similar to PSD but opposite.

### Indefinite

Given a real number $x$, if we say $x$ is not non-negative, then $x$ is negative. However, things are different for matrices. If a matrix is not NND, we can't safely say it's NPD. There exists a big "gray area" in the middle for matrices that are neither NND nor NPD called `indefinite` matrices. For those, $\boldsymbol{x}'\boldsymbol{Ax}$ can be negative and positive depending on $\boldsymbol{x}$.

{{< figure src="matrix_types.png" >}}
