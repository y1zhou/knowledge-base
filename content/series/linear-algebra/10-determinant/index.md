---
title: "Determinant"
date: 2020-10-26T13:25:26-04:00
summary: "The determinant is a very important concept for square matrices, and its properties are key to various other notions such as block matrices and matrix inverses." # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix", "Operation"] # list of related tags

slug: "linear-algebra-determinant"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 100 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "A box of keys." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "zFy6fOPZEu0" # Unsplash ID of the picture
---

The `determinant` is a concept about a square matrix $\boldsymbol{A}$, denoted $|\boldsymbol{A}|$ or $\det(\boldsymbol{A})$. It's a scalar value that can be computed from the matrix.

For a $2 \times 2$ matrix, the determinant may be defined as

$$
\boldsymbol{A} = \begin{pmatrix}
    a & b \\\\ c & d
\end{pmatrix}, \quad |\boldsymbol{A}| = ad - bc
$$

For matrices of higher dimensions, we'll need to first introduce some other concepts.

## Definitions

Let $\boldsymbol{A}$ be an $n \times n$ matrix. The `submatrix` $\boldsymbol{A}_{-i, -j}$ is an $(n-1) \times (n-1)$ matrix by deleting the $i$-th row and $j$-th column from $\boldsymbol{A}$. For example,

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 2 & 3 \\\\
    4 & 5 & 6 \\\\
    7 & 8 & 9
\end{pmatrix}, \quad \boldsymbol{A}\_{-1, -1} = \begin{pmatrix}
    5 & 6 \\\\
    8 & 9
\end{pmatrix}, \quad \boldsymbol{A}\_{-2, -3} = \begin{pmatrix}
    1 & 2 \\\\ 7 & 8
\end{pmatrix}
$$

We can define a submatrix for each $a_{ij}$. We can also find the `co-factor` for each $a_{ij}$:

$$
A_{ij} = (-1)^{i+j} \det(\boldsymbol{A}_{-i, -j})
$$

For example, the co-factor of $a_{11}$ is

$$
A_{11} = (-1)^2 \det \begin{pmatrix}
    5 & 6 \\\\
    8 & 9
\end{pmatrix} = 45 - 48 = -3
$$

Now we can formally define the `determinant` to be

$$
|\boldsymbol{A}| = \det(\boldsymbol{A}) = \sum_{i=1}^n a_{ij}A_{ij},
$$

where $a_{ij}$ is an element of matrix $\boldsymbol{A}$, and $A_{ij}$ is the co-factor of $a_{ij}$. The determinant is a function of $j$ for some fixed $j$, which means we can pick any column of $\boldsymbol{A}$ and the result would be the same. Usually we use the column with the most zeros to simplify the calculation.

Alternatively, we may use any row of $\boldsymbol{A}$ to calculate the determinant:

$$
|\boldsymbol{A}| = \sum_{j=1}^n a_{ij}A_{ij} \quad \text{for any fixed }i
$$

### Calculation of determinant

Suppose

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 1 & 3 \\\\
    1 & 0 & 2 \\\\
    1 & -1 & 1
\end{pmatrix}
$$

To find its determinant, we may use the second row since it contains a zero:

$$
\begin{aligned}
    |\boldsymbol{A}| &= a_{21}A_{21} + a_{22}A_{22} + a_{23}A_{23} \\\\
    &= 1 \times (-1)^3 \det \begin{pmatrix}
        1 & 3 \\\\
        -1 & 1
    \end{pmatrix} + 2 \times (-1)^5 \det \begin{pmatrix}
        1 & 1 \\\\
        1 & -1
    \end{pmatrix} \\\\
    &= -(1 + 3) - 2(-1 - 1) \\\\
    &= -4 + 4 = 0
\end{aligned}
$$

## Properties of the determinant

Knowing how to calculate the determinant is far less important than knowing its properties. For larger matrices we almost always use software to find the deterinant.

1. $|\boldsymbol{A}| = |\boldsymbol{A}'|$.
2. If two rows (columns) are identical, then the determinant is zero.
3. If $\boldsymbol{A}$ is non-singular, it's equivalent to $|\boldsymbol{A}| \neq 0$.
4. If two rows (columns) are exchanged, the sign of the determinant changes, i.e. $|\boldsymbol{A}| \rightarrow -|\boldsymbol{A}|$.
5. $|\boldsymbol{AB}| = |\boldsymbol{A}| |\boldsymbol{B}|$.
6. Because of (5), we also have $|\boldsymbol{AB}| = |\boldsymbol{BA}|$.
7. $|\boldsymbol{AA}'| = |\boldsymbol{A}| |\boldsymbol{A}'| = |\boldsymbol{A}|^2 \geq 0$. Following the same logic we have $|\boldsymbol{A}'\boldsymbol{A}| \geq 0$.
8. If a row (column) is $\boldsymbol{0}$, then $|\boldsymbol{A}| = 0$.
9. If $\boldsymbol{A}$ is diagonal, $\boldsymbol{A} = diag(d_1, \cdots, d_n)$, then

    $$
    |\boldsymbol{A}| = \prod_{i=1}^n d_i = d_1 \times d_2 \times \cdots \times d_n
    $$

10. $|c\boldsymbol{A}| = c^n |\boldsymbol{A}|$.

    $$
    |c\boldsymbol{A}| = |c\boldsymbol{I}_n \boldsymbol{A}| = |c\boldsymbol{I}_n||\boldsymbol{A}| = c^n |\boldsymbol{A}|
    $$

11. If $\boldsymbol{A}$ is triangular, $|\boldsymbol{A}| = \prod_{i=1}^n a_{ii}$. As an example,

    $$
    \left| \begin{pmatrix}
        1 & 2 & 3 \\\\
        0 & 9 & -10 \\\\
        0 & 0 & -3
    \end{pmatrix} \right| = \left| \begin{pmatrix}
        1 & 0 & 0 \\\\
        0 & 9 & 0 \\\\
        0 & 0 & -3
    \end{pmatrix} \right|
    $$

12. Suppose square matrix $\boldsymbol{A}$ is partitioned into

    $$
    \boldsymbol{A} = \begin{bmatrix}
        \boldsymbol{A}\_{11} & \boldsymbol{A}\_{12} \\\\
        \boldsymbol{A}\_{21} & \boldsymbol{A}\_{22}
    \end{bmatrix}
    $$

    such that $\boldsymbol{A}\_{11}$ and $\boldsymbol{A}\_{22}$ are square. If $\boldsymbol{A}\_{12} = 0$ and $\boldsymbol{A}\_{21} = 0$, then $|\boldsymbol{A}| = |\boldsymbol{A}\_{11}| |\boldsymbol{A}\_{22}|$.

    This structure is called `block diagonals` and can obviously be extended to more than two blocks.

## Calculating matrix inverse

To find the inverse of a matrix, we need a final piece of definition called the `adjugate`. The adjugate of a matrix $\boldsymbol{A}$ is the transpose of an $n \times n$ matrix with co-factors

$$
\boldsymbol{A}^* = \\{A_{ij}\\}'
$$

For example,

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 2 & 3 \\\\
    2 & 4 & 5 \\\\
    0 & 1 & 6
\end{pmatrix}, \quad
\\{A_{ij}\\} = \begin{pmatrix}
    19 & -12 & 2 \\\\
    -9 & 6 & -1 \\\\
    -2 & 1 & 0
\end{pmatrix}, \quad
\boldsymbol{A}^* = \begin{pmatrix}
    19 & -9 & -2 \\\\
    -12 & 6 & 1 \\\\
    2 & -1 & 0
\end{pmatrix}
$$

Then the inverse of $\boldsymbol{A}$ can be found by

$$
\boldsymbol{A}^{-1} = \frac{1}{|\boldsymbol{A}|}\boldsymbol{A}^*
$$

Using the same example above, we have

$$
\begin{aligned}
    |\boldsymbol{A}| &= 0 \times (-1)^{3+1} \left| \begin{pmatrix}
        2 & 3 \\\\ 4 & 5
    \end{pmatrix} \right| +
    1 \times (-1)^{3+2} \left| \begin{pmatrix}
        1 & 3 \\\\ 2 & 5
    \end{pmatrix} \right| +
    6 \times (-1)^{3+3} \left| \begin{pmatrix}
        1 & 2 \\\\ 2 & 4
    \end{pmatrix} \right| \\\\
    &= -1 \times (-1) + 6 \times 0 \\\\
    &= 1
\end{aligned}
$$

So the inverse of $\boldsymbol{A}$ is just $\boldsymbol{A}^*$ in this case.

### Properties of matrix inverse

1. $|\boldsymbol{A}^{-1}| = \frac{1}{|\boldsymbol{A}|}$.

    Since $\boldsymbol{AA}^{-1} = \boldsymbol{I}$,

    $$
    |\boldsymbol{AA}^{-1}| = |\boldsymbol{A}| |\boldsymbol{A}^{-1}| = 1
    $$

    Since $|\boldsymbol{A}| \neq 0$, we have

    $$
    |\boldsymbol{A}^{-1}| = \frac{1}{|\boldsymbol{A}|} = |\boldsymbol{A}|^{-1}
    $$

2. If $\boldsymbol{A}$ is non-singular and symmetric, $\boldsymbol{A}^{-1}$ is also non-singular and symmetric, i.e. $(\boldsymbol{A}^{-1})' = \boldsymbol{A}^{-1}$.

    Again since $\boldsymbol{AA}^{-1} = \boldsymbol{I}$, if we take the transpose of both sides:

    $$
    \left(\boldsymbol{AA}^{-1} \right)' = \boldsymbol{I} = (\boldsymbol{A}^{-1})'\boldsymbol{A}' = (\boldsymbol{A}^{-1})'\boldsymbol{A}
    $$

    Thus we have $(\boldsymbol{A}^{-1})' = \boldsymbol{A}^{-1}$. As the inverse matrix is unique, it must be the same inverse that we started with.

3. $(a\boldsymbol{I}_n + b \boldsymbol{J}_n)^{-1} = \frac{1}{a} \left(\boldsymbol{I}_n - \frac{b}{a + nb}\boldsymbol{J}_n \right)$. If $a=0$, the matrix is singular and the inverse doesn't exist. The proof for this is trivial, but it's a very important equation that's frequently used in statistical analyses such as time series analysis.

    Note that the inverse also has a structure of $c\boldsymbol{I} + d\boldsymbol{J}$.

4. $(\boldsymbol{A} + \boldsymbol{CBD})^{-1} = \boldsymbol{A}^{-1} - \boldsymbol{A}^{-1}\boldsymbol{C}(\boldsymbol{B}^{-1} + \boldsymbol{DA}^{-1}\boldsymbol{C})^{-1}\boldsymbol{DA}^{-1}$, assuming all inverses exist. Some special cases for this are:

    1. $(\boldsymbol{A} + \boldsymbol{B})^{-1} = \boldsymbol{A}^{-1} - \boldsymbol{A}^{-1}(\boldsymbol{B}^{-1} + \boldsymbol{A}^{-1})^{-1}\boldsymbol{A}^{-1}$. This is also equal to $\boldsymbol{B}^{-1} - \boldsymbol{B}^{-1}(\boldsymbol{A}^{-1} + \boldsymbol{B}^{-1})^{-1}\boldsymbol{B}^{-1}$.
    2. If $\boldsymbol{B} = \boldsymbol{I}$ and $\boldsymbol{C} = \boldsymbol{c}$, $\boldsymbol{D} = \boldsymbol{d}'$, then

        $$
        \begin{aligned}
            (\boldsymbol{A} + \boldsymbol{cd}')^{-1} &= \boldsymbol{A}^{-1} - \boldsymbol{A}^{-1}\boldsymbol{c}(\boldsymbol{I} + \boldsymbol{d}'\boldsymbol{A}^{-1}\boldsymbol{c})^{-1}\boldsymbol{d}'\boldsymbol{A}^{-1} \\\\
            &= \boldsymbol{A}^{-1} - \frac{\boldsymbol{A}^{-1}\boldsymbol{cd}'\boldsymbol{A}^{-1}}{1 + \boldsymbol{d}'\boldsymbol{A}^{-1}\boldsymbol{c}}
        \end{aligned}
        $$

        To understand this, see that $\boldsymbol{B}$ must be a scalar as its dimensions are limited by $\boldsymbol{c}: n \times 1$ and $\boldsymbol{d}': 1 \times n$.

5. For block-diagonals,

    $$
    \boldsymbol{A} = \begin{pmatrix}
        \boldsymbol{A}\_{11} & \boldsymbol{0} \\\\
        \boldsymbol{0} & \boldsymbol{A}\_{22}
    \end{pmatrix}, \quad
    \boldsymbol{A}^{-1} = \begin{pmatrix}
        \boldsymbol{A}\_{11}^{-1} & \boldsymbol{0} \\\\
        \boldsymbol{0} & \boldsymbol{A}\_{22}^{-1}
    \end{pmatrix}
    $$

6. Another fact for block matrices that's really important in regression analysis (see [graphical models](https://en.wikipedia.org/wiki/Graphical_model) for more information):

$$
\boldsymbol{A} = \left(\begin{array}{c|c}
    \boldsymbol{A}\_{11} & \boldsymbol{A}\_{12} \\\\
    \hline
    \boldsymbol{A}\_{21} & \boldsymbol{A}\_{22}
\end{array}\right)
$$

Let $\boldsymbol{A}\_{11 \cdot 2} = \boldsymbol{A}\_{11} - \boldsymbol{A}\_{12}\boldsymbol{A}\_{22}^{-1}\boldsymbol{A}\_{21}$, and $\boldsymbol{A}\_{22 \cdot 1} = \boldsymbol{A}\_{22} - \boldsymbol{A}\_{21} \boldsymbol{A}\_{11}^{-1} \boldsymbol{A}\_{12}$. The structure of these seemingly complex formulae is important for understanding regression.

$$
\boldsymbol{A}^{-1} = \begin{pmatrix}
    \boldsymbol{A}\_{11 \cdot 2}^{-1} & -\boldsymbol{A}\_{11}^{-1} \boldsymbol{A}\_{12} \boldsymbol{A}\_{22 \cdot 1}^{-1} \\\\
    -\boldsymbol{A}\_{22 \cdot 1}^{-1} \boldsymbol{A}\_{21} \boldsymbol{A}\_{11 \cdot 2} & \boldsymbol{A}\_{22 \cdot 1}^{-1}
\end{pmatrix}
$$
