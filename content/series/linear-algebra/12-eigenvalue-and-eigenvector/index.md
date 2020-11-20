---
title: "Eigenvalues and Eigenvectors"
date: 2020-11-18T18:45:52-05:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Statistics"] # list of related tags

slug: "linear-algebra-eigenvalues-and-eigenvectors"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: true
draft: false

weight: 120 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Dimensions." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "4k3jBXHMEwo" # Unsplash ID of the picture
---

If you're a statistics major (or any science major really), you must have heard of these two words before. In this chapter, we're going to learn about how to find the eigenvalues and eigenvectors of a matrix (eigendecomposition) and singular value decomposition.

## Definitions

A scalar $\lambda$ is an `eigenvalue` of an $n \times n$ matrix $\boldsymbol{A}$ if there exists $\boldsymbol{x} \neq \boldsymbol{0}$ such that

$$
\begin{equation}\label{eq:ev}
    \boldsymbol{Ax} = \lambda \boldsymbol{x}
\end{equation}
$$

where $\boldsymbol{x}$ is called the `eigenvector` corresponding to $\lambda$.

Taking a closer look at Eq.$\eqref{eq:ev}$, we can see that on the LHS we're multiplying matrix $\boldsymbol{A}$ with a vector, and on the RHS we have the scalar multiple of the same vector. If we move the RHS to the left,

$$
\underbrace{(\boldsymbol{A} - \lambda\boldsymbol{I})}\_{n \times n} \boldsymbol{x} = \boldsymbol{0}
$$

The equation implies that $\boldsymbol{x} \in \mathcal{N}(\boldsymbol{A} - \lambda\boldsymbol{I})$. Thus, $\lambda$ is an eigenvalue of $\boldsymbol{A}$ is equivalent to the statement that $\boldsymbol{A} - \lambda\boldsymbol{I}$ is singular, i.e.

$$
\begin{equation}\label{eq:find-ev}
    \det(\boldsymbol{A} - \lambda\boldsymbol{I}) = 0
\end{equation}
$$

Eq.$\eqref{eq:find-ev}$ is generally the equation we solve to find the eigenvalues.

1. If we multiply $\boldsymbol{x}'$ to both sides of Eq.$\eqref{eq:ev}$,

    $$
    \boldsymbol{x}'\boldsymbol{Ax} = \boldsymbol{x}'\lambda\boldsymbol{x} = \lambda \boldsymbol{x}'\boldsymbol{x} \Longrightarrow \lambda = \frac{\boldsymbol{x}'\boldsymbol{Ax}}{\boldsymbol{x}'\boldsymbol{x}}
    $$

    This is the relationship between the eigenvalue and eigenvector.

2. If we multiply a non-zero scalar $c$ to Eq.$\eqref{eq:ev}$,

    $$
    \begin{gathered}
        c\boldsymbol{Ax} = c\lambda\boldsymbol{x} \\\\
        \boldsymbol{A}(c\boldsymbol{x}) = \lambda(c\boldsymbol{x})
    \end{gathered}
    $$

    This means if $\boldsymbol{x}$ is an eigenvector of $\lambda$, then so is $c\boldsymbol{x}$. {{<hl>}}Eigenvectors are not unique!{{</hl>}} It is customary to normalize $\boldsymbol{x}$ such that $\\|\boldsymbol{x}\\| = 1$[^eigenvector-uniqueness].

[^eigenvector-uniqueness]: Note that this still doesn't make eigenvectors unique. For example, if $(1, 0, 0)$ is an eigenvector, then so is $(-1, 0, 0)$. Signs are still ambiguous.

## Calculation

The classical method is to first find the eigenvalues, and then calculate the eigenvectors for each eigenvalue. To find $\lambda$, see that in general Eq.$\eqref{eq:find-ev}$ is an $n^{th}$-degree polynomial of the variable $\lambda$:

$$
a_n\lambda^n + a_{n-1}\lambda^{n-1} + \cdots + a_1\lambda + a_0 = 0
$$

For example, the matrix

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 1 \\\\
    0 & 1
\end{pmatrix}, \quad
\lambda\boldsymbol{I} = \begin{pmatrix}
    \lambda & 0 \\\\
    0 & \lambda
\end{pmatrix}
$$

We have

$$
|\boldsymbol{A} - \lambda\boldsymbol{I}| = \left|\begin{pmatrix}
    1-\lambda & 1 \\\\
    0 & 1-\lambda
\end{pmatrix}\right| = (1-\lambda)^2 = 0
$$

The solution is obviously $\lambda = 1$. Since this comes from a second-order polynomial, we say it has algebraic multiplicity 2. `Algebraic multiplicity` is the number of times $\lambda$ appears as a root of the polynomial.

To find the eigenvector $\boldsymbol{x}$, recall that it's in the nullspace of

$$
\boldsymbol{A} - \lambda\boldsymbol{I} = \begin{pmatrix}
    0 & 1 \\\\
    0 & 0
\end{pmatrix}
$$

So we have

$$
\begin{pmatrix}
    0 & 1 \\\\
    0 & 0
\end{pmatrix}
\begin{pmatrix}
    x_1 \\\\ x_2
\end{pmatrix} = \boldsymbol{0} \Rightarrow x_2 = 0
$$

So $\boldsymbol{x} = (1, 0)'$ is our eigenvector, where $x_1 = 1$ because we want to make the norm 1. Sometimes we would find more than one (say, $k$) eigenvector associated with the same eigenvalue, then we would say the eigenvalue has geometric multiplicity of $k$. The `geometric multiplicity` is the dimension of the eigenspace[^eigenspace].

[^eigenspace]: The eigenspace is the nullspace of $\boldsymbol{A} - \lambda\boldsymbol{I}$. It's the space spanned by the eigenvectors of $\boldsymbol{A}$ for $\lambda$.

In general, the algebraic multiplicity and geometric multiplicity of an eigenvalue can differ, but the geometric multiplicity can **never exceed** the algebraic multiplicity.

### Notes

1. **Diagonal matrices**: it's easy to calculate the eigenvalues of diagonal matrices. For example,

    $$
    \begin{gathered}
        \boldsymbol{D} = diag(3, 1, -5, -5, 0) \\\\
        \boldsymbol{D} - \lambda\boldsymbol{I} = diag(3-\lambda, 1-\lambda, -5-\lambda, -5-\lambda, -\lambda) \\\\
        |\boldsymbol{D} - \lambda\boldsymbol{I}| = \prod_{i=1}^5 d_i = -(3-\lambda)(1-\lambda)(-5-\lambda)^2\lambda = 0
    \end{gathered}
    $$

    The solutions for $\lambda$ are 0, 1, 3, -5 where -5 has algebraic multiplicity 2. Often the eigenvalues of an diagonal matrix are the diagonals themselves.

2. **Transpose**: eigenvalues of $\boldsymbol{A}'$ are the same as eigenvalues of $\boldsymbol{A}$. $|\boldsymbol{A}' - \lambda\boldsymbol{I}| = 0$ if $\lambda$ is an eigenvalue of $\boldsymbol{A}$ because $|\boldsymbol{A}' - \lambda\boldsymbol{I}| = |\boldsymbol{A} - \lambda\boldsymbol{I}|$.

3. **Powers**: for powers of $\boldsymbol{A}$, $\boldsymbol{A}^k$, we have from Eq.$\eqref{eq:ev}$

    $$
    \boldsymbol{A}^2\boldsymbol{x} = \boldsymbol{A}\lambda\boldsymbol{x} = \lambda\boldsymbol{Ax} = \lambda^2\boldsymbol{x},
    $$

    which means $(\lambda^2, \boldsymbol{x})$ is a pair of eigenvalue and eigenvector for $\boldsymbol{A}^2$. Doing this recursively and we'll get $(\lambda^k, \boldsymbol{x})$ is an eigenvalue and eigenvector pair for $\boldsymbol{A}^k$.

4. **Inverse**: it's natural to pre-multiply $\boldsymbol{A}^{-1}$ to both sides of Eq.$\eqref{eq:ev}$ and get

    $$
    \boldsymbol{x} = \lambda\boldsymbol{A}^{-1}\boldsymbol{x} \Rightarrow \boldsymbol{A}^{-1}\boldsymbol{x} = \frac{1}{\lambda}\boldsymbol{x}
    $$

    assuming $\lambda \neq 0$. Thus $(\frac{1}{\lambda}, \boldsymbol{x})$ is an eigenvalue and eigenvector pair for $\boldsymbol{A}^{-1}$.

    In fact, all eigenvalues of a non-singular matrix are non-zero. Because if $\lambda = 0$, then from Eq.$\eqref{eq:ev}$ we have $\boldsymbol{Ax} = \boldsymbol{0}$, meaning that $\boldsymbol{x} \in \mathcal{N}(\boldsymbol{A})$ and thus the dimension of the nullspace of $\boldsymbol{A}$ is at least one, which contradicts the fact that $\boldsymbol{A}$ is non-singular.

### Examples

-   In the first example, let

    $$
    \boldsymbol{A} = \begin{pmatrix}
        a & b \\\\
        c & d
    \end{pmatrix}, \quad
    \boldsymbol{A} - \lambda\boldsymbol{I} = \begin{pmatrix}
        a - \lambda & b \\\\
        c & d - \lambda
    \end{pmatrix}
    $$

    So the determinant is

    $$
    \begin{aligned}
        \det(\boldsymbol{A} - \lambda\boldsymbol{I}) &= (a - \lambda)(d-\lambda) - bc \\\\
        &= \lambda^2 - (a+d)\lambda + ad - bc \\\\
        &= (\lambda - \lambda_1)(\lambda - \lambda_2) = 0
    \end{aligned}
    $$

    We can actually get

    $$
    \begin{gathered}
        \lambda_1 + \lambda_2 = a + d = tr(\boldsymbol{A}) \\\\
        \lambda_1\lambda_2 = ad - bc = \det(\boldsymbol{A})
    \end{gathered}
    $$

    This holds true for not just $2 \times 2$ matrices. This also tells us that $\boldsymbol{A}$ is singular if and only if 0 is an eigenvalue.

-   In the second example, let

    $$
    \boldsymbol{A} = \begin{pmatrix}
        2 & 2 & 0 \\\\
        2 & 1 & 1 \\\\
        -7 & 2 & -3
    \end{pmatrix}
    $$

    To find the eigenvalues,

    $$
    \begin{aligned}
        |\boldsymbol{A} - \lambda\boldsymbol{I}_3| &= \lambda^3 - 13\lambda + 12 \\\\
        &= (\lambda-1)(\lambda+4)(\lambda-3) = 0
    \end{aligned}
    $$

    We have $\lambda_1 = 3$, $\lambda_2 = 1$, and $\lambda_3 = -4$[^eigenvalue-order]. Because there are three unique eigenvalues, we know the dimension of the eigenspace would be 1. Now for the eigenvector(s) associated with $\lambda = 1$,

    [^eigenvalue-order]:
        In practice, we usually put the eigenvalues in decreasing order, so

        $$
        \lambda_1 \geq \lambda_2 \geq \lambda_3 \geq \cdots
        $$

    $$
    \begin{gathered}
        \boldsymbol{Ax} = \boldsymbol{x} \Rightarrow (\boldsymbol{A} - \boldsymbol{I})\boldsymbol{x} = \boldsymbol{0} \\\\
        (\boldsymbol{A} - \boldsymbol{I})\boldsymbol{x} = \begin{pmatrix}
            1 & 2 & 0 \\\\
            2 & 0 & 1 \\\\
            -7 & 2 & -4
        \end{pmatrix}
        \begin{pmatrix}
            x_1 \\\\ x_2 \\\\ x_3
        \end{pmatrix} \\\\
        \begin{cases}
            x_1 + 2x_2 = 0 \\\\
            2x_1 + x_3 = 0 \\\\
            -7x_1 + 2x_2 - 4x_3 = 0
        \end{cases} \Rightarrow
        \begin{cases}
            x_1 = -2x_2 \\\\
            -4x_2 + x_3 = 0
        \end{cases}
    \end{gathered}
    $$

    If $x_2 = 1$, we have $\boldsymbol{x} = (-2, 1, 4)'$, which can be normalized to $\boldsymbol{x} = \frac{1}{\sqrt{21}}(-2, 1, 4)'$. This is the eigenvector of $\boldsymbol{A}$ that corresponds to $\lambda = 1$.

    Repeating the same procedure and we can find the corresponding eigenvectors for $\lambda = -4$ to be $\boldsymbol{x} = \frac{1}{3}(-2, -1, 2)'$ and for $\lambda = 3$ to be $\boldsymbol{x} = \frac{1}{\sqrt{179}}(1, -3, 13)'$. We can see that the three vectors are linearly independent.

-   The third example is an symmetric matrix:

    $$
    \boldsymbol{A} = \begin{pmatrix}
        1 & 2 & 2 \\\\
        2 & 1 & 2 \\\\
        2 & 2 & 1
    \end{pmatrix}
    $$

    First we find the eigenvalues:

    $$
    \begin{gathered}
        |\boldsymbol{A} - \lambda\boldsymbol{I}| = (\lambda+1)(\lambda-5) = 0 \\\\
        \lambda_1 = 5, \lambda_2 = -1
    \end{gathered}
    $$

    where $\lambda_2 = -1$ has algebraic multiplicity 2. Its corresponding eigenvector can be found with

    $$
    (\boldsymbol{A} - \lambda_2\boldsymbol{I})\boldsymbol{x} = 2\boldsymbol{J}_3\boldsymbol{I} = \boldsymbol{0}
    $$

    and we get $x_1 + x_2 + x_3 = 0$. The nullspace of $\boldsymbol{A} - \lambda_1\boldsymbol{I}$ is $\\{\boldsymbol{x} \mid x_1 + x_2 + x_3 = 0 \\}$ and its dimension is 2. Any vector in this subspace is an eigenvector for $\lambda_1$, and this space is called the eigenspace.

    It is customary to provide an orthonormal basis of the eigenspace. In this case,

    $$
    \boldsymbol{x}_1 = \frac{1}{\sqrt{6}}(2, -1, -1)', \boldsymbol{x}_2 = \frac{1}{\sqrt{2}}(0, 1, -1)'
    $$

    Similarly for $\lambda_1 = 5$ we can find $\boldsymbol{x}_3 = \frac{1}{\sqrt{3}}(1, 1, 1)'$.
