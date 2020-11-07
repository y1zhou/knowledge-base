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

## Properties

1.  If $\boldsymbol{A}$ is NND/PD/PSD, then $-\boldsymbol{A}$ is NPD/ND/NSD, respectively.
2.  If $\boldsymbol{D}$ is a diagonal matrix $\boldsymbol{D} = diag(d_1, \cdots, d_n)$, then

    $$
    d_i \geq 0 \\, \forall i \Longleftrightarrow \boldsymbol{D} \text{ is NND}
    $$

    This is because

    $$
    \boldsymbol{x}'\boldsymbol{Dx} = d_1x_1^2 + d_2x_2^2 + \cdots + d_nx_n^2 \geq 0
    $$

    Following the same logic, we have

    $$
    d_i > 0 \\, \forall i \Longleftrightarrow \boldsymbol{D} \text{ is PD}
    $$

3.  Multiplying a positive scalar to $\boldsymbol{A}$ does not change the definiteness.
4.  If $\boldsymbol{A}$ and $\boldsymbol{B}$ are both NND, then $\boldsymbol{A} + \boldsymbol{B}$ is NND[^nnd-other-cases].

    $$
    \boldsymbol{x}'(\boldsymbol{A} + \boldsymbol{B})\boldsymbol{x} = \boldsymbol{x}'\boldsymbol{Ax} + \boldsymbol{x}'\boldsymbol{Bx} \geq 0
    $$

    Similarly if $\boldsymbol{A}$ and $\boldsymbol{B}$ are both PD, then $\boldsymbol{A} + \boldsymbol{B}$ is PD.

    However, if $\boldsymbol{A}$ and $\boldsymbol{B}$ are both PSD, $\boldsymbol{A} + \boldsymbol{B}$ is not necessarily PSD. $\boldsymbol{A} + \boldsymbol{B}$ is NND, but it could be PD or PSD.

    Finally, if $\boldsymbol{A}$ is PSD and $\boldsymbol{B}$ is PD, then $\boldsymbol{A} + \boldsymbol{B}$ is PD.

[^nnd-other-cases]: These conclusions cannot be generalized to $\boldsymbol{A} - \boldsymbol{B}$ or $\boldsymbol{AB}$.

5. {{<hl>}}A PD matrix is non-singular.{{</hl>}} An important corollary of this is that a PSD matrix is singular. This is one of the reasons why we put an emphasis on differentiating these two cases -- singularity comes with a lot of properties!

    To prove this, suppose $\boldsymbol{A}$ is singular. This means that the nullspace of $\boldsymbol{A}$ is not empty, i.e. there exists $\boldsymbol{x} \neq \boldsymbol{0}$ such that $\boldsymbol{x} \in \mathcal{N}(\boldsymbol{A})$. We then have

    $$
    \boldsymbol{Ax} = \boldsymbol{0} \Rightarrow \boldsymbol{x}'\boldsymbol{Ax} = 0
    $$

    We have found a non-zero $\boldsymbol{x}$ that can make the quadratic form zero, which is contradictory to "$\boldsymbol{A}$ is PD".

6. Suppose $\boldsymbol{A}: n \times n$ and $\boldsymbol{P}: n \times m$.

    1. If $\boldsymbol{A}$ is NND, then $\boldsymbol{P}'\boldsymbol{AP}$ is NND. This can be easily proven using the quadratic form:

        $$
        \boldsymbol{x}'\boldsymbol{P}'\boldsymbol{APx} = (\boldsymbol{Px})'\boldsymbol{A}(\boldsymbol{Px}) \Rightarrow \boldsymbol{y}'\boldsymbol{Ay} \geq 0
        $$

    2. If $\boldsymbol{A}$ is PD and $r(\boldsymbol{P}) = m$, then $\boldsymbol{P}'\boldsymbol{AP}$ is PD.

        To show this, consider when $\boldsymbol{x}'\boldsymbol{P}'\boldsymbol{APx} = 0$. The fact that $\boldsymbol{A}$ is PD leads to $\boldsymbol{Px} = \boldsymbol{0}$. Since $r(\boldsymbol{P}) = m$, the $m \times m$ matrix $\boldsymbol{P}'\boldsymbol{P}$ is non-singular, and $\boldsymbol{P}'\boldsymbol{Px} = \boldsymbol{0}_m$ leads to $\boldsymbol{x} = \boldsymbol{0}$.

        Similarly we can show that if $\boldsymbol{A}$ is PSD and $r(\boldsymbol{P}) = m$, then $\boldsymbol{P}'\boldsymbol{AP}$ is PSD.

    3. If $\boldsymbol{A}$ is NND and $r(\boldsymbol{P}) < m$, then $\boldsymbol{P}'\boldsymbol{AP}$ is PSD. From the properties of matrix rank we know that $r(\boldsymbol{P}'\boldsymbol{AP}) < m$, so it's singular and can't be PD.
