---
title: "Matrices"
date: 2020-08-26T15:14:34-04:00
summary: "Matrix algebra plays an important role in many areas of statistics, such as linear statistical models and multivariate analysis. In this chapter we introduce basic terminology and some basic matrix operations. We also introduce some basic types of matrices." # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix", "Operations"] # list of related tags

slug: "linear-algebra-matrices"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 10 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Matrix code." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "iar-afB0QQw" # Unsplash ID of the picture
---

### Basic terminology

A `matrix` is a rectangular array of real numbers. It takes the general form

$$
\boldsymbol{A}\_{m \times n} = \begin{pmatrix}
    a_{11} & a_{12} & \cdots & a_{1n} \\\\
    a_{21} & a_{22} & \cdots & a_{2n} \\\\
    \vdots & \vdots & & \vdots \\\\
    a_{m1} & a_{m2} & \cdots & a_{mn}
\end{pmatrix}
$$

The convention is to represent a matrix with boldface uppercase letters, and the `dimensions` (size) of the matrix $m$ by $n$ is always rows by columns.

If $m=1$, $A_{1 \times n} = [a_{11}, \cdots, a_{1n}]$ is a `row vector` denoted by $\underline{a}^\prime$. Similarly if $n=1$, we have a `column vector` that's denoted by $\underline{a}$.

If two matrices $\boldsymbol{A} = \\{a_{ij}\\}$ and $\boldsymbol{B} = \\{b_{ij}\\}$ have the same dimensions, the following basic operations apply:

$$
\begin{gathered}
    \boldsymbol{A} + \boldsymbol{B} = \\{a_{ij} + b_{ij}\\} \\\\
    \boldsymbol{A} - \boldsymbol{B} = \\{a_{ij} - b_{ij}\\} \\\\
    c\boldsymbol{A} = \\{ca_{ij}\\}
\end{gathered}
$$

For `matrix multiplication` (product), $\boldsymbol{AB}$ is defined only when the number of columns of $\boldsymbol{A}$ and the number of rows of $\boldsymbol{B}$ are equal. If $\boldsymbol{A}: m \times n$ and $\boldsymbol{B}: n \times p$, the size of the product $\boldsymbol{C} = \boldsymbol{AB}$ is $m \times p$. The $(i, j)^{th}$ element of $\boldsymbol{C}$ is given by

$$
c_{ij} = \sum_{k=1}^n a_{ik} b_{kj} = a_{i1}b_{1j} + a_{i2}b_{2j} + \cdots + a_{in}b_{nj}
$$

which is the inner product of the $i^{th}$ row of $\boldsymbol{A}$ and the $j^{th}$ column of $\boldsymbol{B}$.

In general, matrix multiplication is not commutative, that is $\boldsymbol{AB}$ is not necessarily identical to $\boldsymbol{BA}$. We refer to the matrix product $\boldsymbol{AB}$ as $\boldsymbol{A}$ is `premultiplied` to $\boldsymbol{B}$, or $\boldsymbol{B}$ is `postmultiplied` to $\boldsymbol{A}$.

Matrix multiplication is however associative, meaning that $\boldsymbol{A}(\boldsymbol{BC}) = (\boldsymbol{AB})\boldsymbol{C}$, assuming that $\boldsymbol{A}: m \times n$, $\boldsymbol{B}: n \times p$ and $\boldsymbol{C}: p \times q$.

One way to prove this is to take the $(i, j)^{th}$ element on the left hand side is equal to the $(i, j)^{th}$ element on the right hand side. Another way is to take the left hand side, go through some transformations and get the right hand side. There's several other ways, such as showing $\boldsymbol{A} - \boldsymbol{B} = \boldsymbol{0}$, or $(\boldsymbol{A} - \boldsymbol{B})^\prime(\boldsymbol{A} - \boldsymbol{B}) = \boldsymbol{0}$.

**Proof:** The $(i, j)^{th}$ element on the LHS is

$$
\begin{aligned}
    \sum_{k=1}^p (\boldsymbol{AB})_{ik} c_{kj} &= \sum_{k=1}^p \left( \sum_{l=1}^n a_{il} b_{lk} \right) c_{kj} \\\\
    &= \sum_{l=1}^n \sum_{k=1}^p a_{il} b_{lk} c_{kj} \\\\
    &= \sum_{l=1}^n a_{il} \sum_{k=1}^p b_{lk} c_{kj} \\\\
    &= \boldsymbol{A}(\boldsymbol{BC})_{ij}
\end{aligned}
$$

Matrix multiplication is distributive, that is $\boldsymbol{A}(\boldsymbol{B} + \boldsymbol{C}) = \boldsymbol{AB} + \boldsymbol{AC}$. This property can be proved in a similar fashion.

### Transpose

The `transpose` of an $m \times n$ matrix $\boldsymbol{A}$, denoted $\boldsymbol{A}^\prime$ or $\boldsymbol{A}^T$, is the $n \times m$ matrix whose $(i, j)^{th}$ element is the $(j, i)^{th}$ element of $\boldsymbol{A}$. The columns of $\boldsymbol{A}$ becomes the rows of $\boldsymbol{A}^\prime$.

For any matrix $\boldsymbol{A}$, $(\boldsymbol{A}^\prime)^\prime = \boldsymbol{A}$.

For any two matrices $\boldsymbol{A}$ and $\boldsymbol{B}$, $(\boldsymbol{AB})^\prime = \boldsymbol{B}^\prime \boldsymbol{A}^\prime$. In general,

$$
(\boldsymbol{A}_1 \boldsymbol{A}_2 \cdots \boldsymbol{A}_k)^\prime = \boldsymbol{A}_k^\prime \cdots \boldsymbol{A}_1^\prime
$$

**Proof:** Let's say $\boldsymbol{A}: m \times n$ and $\boldsymbol{B}: n \times p$. The $(i, j)^{th}$ element of $(\boldsymbol{AB})^\prime$ is equal to the $(j, i)^{th}$ element of $\boldsymbol{AB}$, which is

$$
\sum_{l=1}^n a_{jl}b_{li} = \sum_{l=1}^n b_{li}a_{jl} = \sum_{l=1}^n \boldsymbol{B}_{il}^\prime \boldsymbol{A}_{lj}^\prime
$$

This is the $(i, j)^{th}$ element of $\boldsymbol{B}^\prime \boldsymbol{A}^\prime$.

### Basic types of matrices

#### Square matrix

A matrix having the same number of rows as columns is called a `square matrix`. For square matrices we can talk about `power matrices`. Let $\boldsymbol{A}$ be an $m \times m$ square matrix. The power matrix is defined as

$$
\boldsymbol{A}^k = \underbrace{\boldsymbol{A} \cdot \boldsymbol{A} \cdots \boldsymbol{A}}_k
$$

If $\boldsymbol{A}: m \times n$, $\boldsymbol{A}^\prime \boldsymbol{A}$ is always an $n \times n$ square matrix, and $\boldsymbol{A}\boldsymbol{A}^\prime$ is always an $m \times m$ square matrix.

#### Symmetric matrices

A matrix $\boldsymbol{A}$ is said to be `symmetric` if $\boldsymbol{A}^\prime = \boldsymbol{A}$. A symmetric matrix has to be square whose $ij$th element equals its $ji$th element.

The `diagonal` of an $n \times n$ square matrix refers to the $n$ elements $a_{11}, a_{22}, \cdots, a_{nn}$. The off-diagonal elements are symmetric along the diagonal for symmetric matrices.

If a matrix is symmetric, we have

$$
\boldsymbol{A}^2 = \boldsymbol{A}\boldsymbol{A} = \boldsymbol{A}^\prime \boldsymbol{A} = \boldsymbol{A}\boldsymbol{A}^\prime = (\boldsymbol{A}^\prime)^2 = (\boldsymbol{A}^2)^\prime
$$

If $\boldsymbol{A}$ is square but not symmetric, $\frac{\boldsymbol{A} + \boldsymbol{A}^\prime}{2}$ is always symmetric with the same diagonal elements of $\boldsymbol{A}$.

#### Diagonal matrices

A `diagonal matrix` is a square matrix whose off-diagonal elements are all equal to 0. It's often denoted

$$
\boldsymbol{D} = diag(d_1, \cdots, d_n) = \begin{pmatrix}
    d_1 & 0 & \cdots & 0 \\\\
    0 & d_2 & \cdots & 0 \\\\
    \vdots & & \ddots & \\\\
    0 & 0 & \cdots & d_n
\end{pmatrix}
$$

The power matrix of a diagonal matrix is $\boldsymbol{D}^k = diag(d_1^k, \cdots, d_n^k)$.

For any $m \times n$ matrix $\boldsymbol{A}$ and diagonal matrix $\boldsymbol{D} = diag(d_1, \cdots, d_m)$,

$$
\boldsymbol{DA} = \begin{pmatrix}
    d_1 & & \\\\
    & \ddots & \\\\
    & & d_m
\end{pmatrix}
\begin{pmatrix}
    a_{11} & \cdots & a_{1n} \\\\
    \vdots & \ddots & \\\\
    a_{m1} & \cdots & a_{mn}
\end{pmatrix} =
\begin{pmatrix}
    d_1 a_{11} & d_1 a_{12} & \cdots & d_1 a_{1n} \\\\
    d_2 a_{21} & d_2 a_{22} & \cdots & d_2 a_{2n} \\\\
    \vdots & \vdots & \ddots & \vdots \\\\
    d_m a_{m1} & d_m a_{m2} & \cdots & d_m a_{mn} \\\\
\end{pmatrix}
$$

The effect of premultiplying $\boldsymbol{A}$ by $\boldsymbol{D}$ is to multiply each element of the $i$th row of $\boldsymbol{A}$ by the $i$th diagonal element of $\boldsymbol{D}$. Similarly, the effect of postmultiplying $\boldsymbol{A}$ by a diagonal matrix $\boldsymbol{D}$ of order $n$ is to multiply each element of the $j$th column of $\boldsymbol{A}$ by the $j$th diagonal element of $\boldsymbol{D}$.

The effect of multiplying an $m \times n$ matrix $\boldsymbol{A}$ by a scalar $k$ is the same as that of pre/post-multiplying $\boldsymbol{A}$ by the $m \times m$ or $n \times n$ matrix $diag(k, \cdots, k)$.

#### Identity matrices

A $n \times n$ diagonal matrix with all diagonal elements equal to 1 and off-diagonal elements equal to 0 is called an `identity matrix`, denoted $\boldsymbol{I}_n$. Clearly for an arbitrary matrix $\boldsymbol{A}$,

$$
\boldsymbol{IA} = \boldsymbol{AI} = \boldsymbol{A}
$$

**Proof:** the $ij$th element of the LHS is

$$
\sum_{k=1}^n a_{ik}(\boldsymbol{I}_n)_{kj} = a_{i1}(\boldsymbol{I}_n)_{1j} + \cdots + a_{in}(\boldsymbol{I}_n)_{nj} = a_{ij}
$$

#### Matrix of ones

The symbol $\boldsymbol{J}_{mn}$ is used to represent an $m \times n$ matrix whose elements are all equal to 1. This matrix is very frequently used and it's important to understand it's properties.

-   If $n=1$, $\boldsymbol{J}_{m1} = \boldsymbol{1}_m$.
-   If $m=1$, $\boldsymbol{J}_{1n} = \boldsymbol{1}_n^\prime$.
-   If $\boldsymbol{A}: m \times n$, $\boldsymbol{A1}_n$ gives the sums of the rows in $\boldsymbol{A}$.
-   $\boldsymbol{J}\_{n}\boldsymbol{J}\_{n} = \boldsymbol{J}\_n^2 = n\boldsymbol{J}_n$. Similarly $\boldsymbol{J}_n^3 = n^2\boldsymbol{J}_n$.
-   $\boldsymbol{J}\_{mn}\boldsymbol{J}\_{np} = n\boldsymbol{J}\_{mp}$.

#### Null matrices

A matrix whose elements are all equal to 0 is called a `null matrix` or zero matrix. A null matrix is denoted by the symbol $\boldsymbol{0}$. Note that $\boldsymbol{AB} = \boldsymbol{0}$ does **not** imply $\boldsymbol{A} = \boldsymbol{0}$ or $\boldsymbol{B} = \boldsymbol{0}$!

#### Triangular matrices

If all the elements of a _square_ matrix that are located below and to the left of the diagonal are 0, the matrix is called an `upper triangular matrix`. The general form is

$$
\begin{pmatrix}
    a_{11} & a_{12} & a_{13} & \cdots & a_{1n} \\\\
    0 & a_{22} & a_{23} & \cdots & a_{2n} \\\\
    0 & 0 & a_{33} & \cdots & a_{3n} \\\\
    \vdots & \vdots & & \ddots & \vdots \\\\
    0 & 0 & 0 & & a_{nn}
\end{pmatrix}
$$

More formally, an $n \times n$ matrix $\boldsymbol{A} = \\{a_{ij}\\}$ is upper triangular if $a_{ij} = 0$ for $j < i$, and is lower triangular if $a_{ij} = 0$ for $j > i$. A square matrix is a diagonal matrix if and only if it's both upper and lower triangular.

**Lemma:** If $\boldsymbol{A}$ and $\boldsymbol{B}$, of the same size, are both upper-triangular, then $\boldsymbol{AB}$ is also upper-triangular. Further, if $\boldsymbol{A}$ and $\boldsymbol{B}$ are either both upper-triangular or both lower-triangular, then the $i$th diagonal element of $\boldsymbol{AB}$ is the product $a_{ii}b_{ii}$ of the $i$th diagonal elements of $\boldsymbol{A}$ and $\boldsymbol{B}$.

#### Row and column vectors

A matrix that has only one column is called a `column vector`, and a matrix that has only one row is called a `row vector`. Column vectors can be denoted:

$$
\underline{a} = \boldsymbol{a} = \begin{pmatrix}
    a_1 \\\\ \vdots \\\\ a_n
\end{pmatrix} = (a_1, \cdots, a_n)^\prime
$$

Row vectors are represented as the transposes of columns vectors, e.g. $\boldsymbol{b}^\prime$. If $\boldsymbol{a}$ and $\boldsymbol{b}$ are both $n \times 1$ column vectors, the `inner product` between $\boldsymbol{a}, \boldsymbol{b}$ is $\boldsymbol{a}' \boldsymbol{b}$:

$$
\boldsymbol{a}'\boldsymbol{b} = a_1b_1 + \cdots + a_n b_n = \sum_{i=1}^n a_ib_i
$$

> The inner product is a general operation, but in this class we only use the inner product of two column vectors.

The `outer product` between $\boldsymbol{a}: m \times 1$ and $\boldsymbol{b}': n \times 1$ is $\boldsymbol{ab}'$:

$$
\boldsymbol{ab}' = \begin{pmatrix}
    a_1b_1 & a_1b_2 & \cdots & a_1b_n \\\\
    a_2b_1 & a_2b_2 & \cdots & a_2b_n \\\\
    \vdots & \vdots & \ddots & \vdots \\\\
    a_mb_1 & a_mb_2 & \cdots & a_mb_n
\end{pmatrix} \triangleq \boldsymbol{C}: m \times n
$$

where $\boldsymbol{C}$ is the resulting matrix and $c_{ij} = a_i b_j$.

If $\boldsymbol{A}: m \times n$ and $\boldsymbol{x}: n \times 1$, we can express $\boldsymbol{A}$ using its row and column vectors:

$$
\boldsymbol{A} = (\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n) = \begin{pmatrix}
    \boldsymbol{\alpha}_1^\prime \\\\ \vdots \\\\ \boldsymbol{\alpha}_m^\prime
\end{pmatrix}
$$

The product $\boldsymbol{Ax}$ is

$$
\begin{aligned}
    \boldsymbol{Ax} &= \begin{pmatrix}
        \boldsymbol{\alpha}_1^\prime \\\\ \vdots \\\\ \boldsymbol{\alpha}_m^\prime
    \end{pmatrix}
    \begin{pmatrix}
        \boldsymbol{x}_1 \\\\ \vdots \\\\ \boldsymbol{x}_n
    \end{pmatrix} =
    \begin{pmatrix}
        \boldsymbol{\alpha}_1^\prime \boldsymbol{x} \\\\ \vdots \\\\ \boldsymbol{\alpha}_m^\prime \boldsymbol{x}
    \end{pmatrix} \\\\
    &= (\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n) \begin{pmatrix}
        \boldsymbol{x}_1 \\\\ \vdots \\\\ \boldsymbol{x}_n
    \end{pmatrix} = \cdots = \sum\_{j=1}^n x_j \boldsymbol{a}_j
\end{aligned}
$$

which is a linear combination of **columns of $\boldsymbol{A}$** with $\boldsymbol{x}$ as coefficients.

Similarly, if $\boldsymbol{y}: m \times 1$,

$$
\boldsymbol{y}'\boldsymbol{A} = \sum_{j=1}^m y_j \boldsymbol{\alpha}_j
$$

which is a linear combination of **rows of $\boldsymbol{A}$** with $\boldsymbol{y}$ as coefficients.
