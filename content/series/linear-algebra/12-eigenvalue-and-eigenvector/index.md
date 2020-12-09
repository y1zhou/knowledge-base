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

## Similar matrices

{{<alert info>}}
This section is not much related to the discussions above, but can come useful in some cases.
{{</alert>}}

If matrix $\boldsymbol{B} = \boldsymbol{C}^{-1}\boldsymbol{AC}$, then we say $\boldsymbol{B}$ and $\boldsymbol{A}$ are `similar`. This is equivalent to

$$
\boldsymbol{CB} = \boldsymbol{AC}
$$

We've actually seen one of these examples before when we were discussing symmetric matrix [diagonalization]({{<relref "../11-quadratic-form/index.md#theorem-1">}}). If $\boldsymbol{A}$ is symmetric, it can be written as $\boldsymbol{A} = \boldsymbol{P}^{-1}\boldsymbol{DP}$. A symmetric matrix is similar to a diagonal matrix.

It is implied that both matrices are square. Some properties of similar matrices are:

1. **Rank**: $r(\boldsymbol{B}) = r(\boldsymbol{A})$.
2. **Trace**: The trace of similar matrices are also equal.

    $$
    tr(\boldsymbol{B}) = tr(\boldsymbol{C}^{-1}\boldsymbol{AC}) = tr(\boldsymbol{ACC}^{-1}) = tr(\boldsymbol{A})
    $$

3. **Determinant**: assuming both matrices are square,

    $$
    |\boldsymbol{B}| = |\boldsymbol{C}^{-1}\boldsymbol{AC}| = |\boldsymbol{C}^{-1}| |\boldsymbol{A}| |\boldsymbol{C}| = |\boldsymbol{A}|
    $$

    Here we used the fact that when $|\boldsymbol{C}| \neq 0$, $|\boldsymbol{C}^{-1}| = \frac{1}{|\boldsymbol{C}|}$.

4. **Eigenvalues**: the eigenvalues of $\boldsymbol{B}$ makes $|\boldsymbol{B} - \lambda\boldsymbol{I}| = 0$.

    $$
    \begin{aligned}
        |\boldsymbol{B} - \lambda\boldsymbol{I}| &= |\boldsymbol{C}^{-1}\boldsymbol{AC} - \lambda\boldsymbol{I}| \\\\
        &= |\boldsymbol{C}^{-1}\boldsymbol{AC} - \lambda\boldsymbol{C}^{-1}\boldsymbol{IC}| \\\\
        &= \left| \boldsymbol{C}^{-1} (\boldsymbol{A} - \lambda\boldsymbol{I})\boldsymbol{C} \right| \\\\
        &= | \boldsymbol{C}^{-1}| |\boldsymbol{A} - \lambda\boldsymbol{I}| |\boldsymbol{C}| \\\\
        &= |\boldsymbol{A} - \lambda\boldsymbol{I}|
    \end{aligned}
    $$

    So $\boldsymbol{A}$ and $\boldsymbol{B}$ share the same set of eigenvalues.

5. **Eigenvectors**: keep in mind that $\boldsymbol{B} = \boldsymbol{C}^{-1}\boldsymbol{AC}$.

    $$
    \begin{gathered}
        \boldsymbol{Ax} = \lambda\boldsymbol{x} \\\\
        \boldsymbol{C}^{-1}\boldsymbol{Ax} = \lambda\boldsymbol{C}^{-1}\boldsymbol{x} \\\\
        \boldsymbol{C}^{-1}\boldsymbol{ACC}^{-1}\boldsymbol{x} = \lambda\boldsymbol{C}^{-1}\boldsymbol{x} \\\\
        \boldsymbol{BC}^{-1}\boldsymbol{x} = \lambda\boldsymbol{C}^{-1}\boldsymbol{x}
    \end{gathered}
    $$

    This means if $\boldsymbol{x}$ is an eigenvector of $\boldsymbol{A}$, then $\boldsymbol{C}^{-1}\boldsymbol{x}$ is an eigenvector of $\boldsymbol{B}$. The eigenvectors are **not** shared but can be easily found.

## Important theorems

We are going to establish two very important theorems about eigenvalues and eigenvectors.

### Theorem 1

{{<hl>}}Eigenvectors that are associated with distinct eigenvalues are linearly independent.{{</hl>}}

**Proof**: let $\boldsymbol{v}_1, \cdots, \boldsymbol{v}_n$ be eigenvectors of $\boldsymbol{A}$ with distinct $\lambda_1, \cdots, \lambda_n$, i.e. $\lambda_i \neq \lambda_j$. Suppose $\boldsymbol{v}_i$'s are linearly dependent, consider

$$
\sum_{i=1}^n \alpha_i\boldsymbol{v}_i = \boldsymbol{0}
$$

There exists at least one $\alpha_i \neq 0$. For convenience, we rearrange the coefficients so that $\alpha_n \neq 0$. This means we can express $\boldsymbol{v}_n$ as a linear combination of the other $\boldsymbol{v}_i$'s.

$$
\begin{equation}\label{eq:vn}
    \boldsymbol{v}_n = b_1 \boldsymbol{v}_1 + \cdots + b\_{n-1}\boldsymbol{v}\_{n-1}, \quad b_i = \frac{\alpha_i}{\alpha_n}
\end{equation}
$$

Because $\boldsymbol{v}_n \neq \boldsymbol{0}$, not all $b_i$'s are zeros. Now we multiply $\boldsymbol{A}$ to Eq.$\eqref{eq:vn}$:

$$
\begin{align}
    \boldsymbol{Av}_n &= b_1 \boldsymbol{Av}_1 + \cdots + b\_{n-1}\boldsymbol{Av}\_{n-1} \\\\
    \lambda_n\boldsymbol{v}_n &= b_1 \lambda_1\boldsymbol{v}_1 + \cdots + b\_{n-1} \lambda\_{n-1}\boldsymbol{v}\_{n-1} \label{eq:avn}
\end{align}
$$

Then we multiply $\lambda_n$ to Eq.$\eqref{eq:vn}$ to get

$$
\begin{equation}\label{eq:lambdavn}
    \lambda_n\boldsymbol{v}_n = b_1\lambda_n\boldsymbol{v}_1 + \cdots + b\_{n-1}\lambda_n\boldsymbol{v}\_{n-1}
\end{equation}
$$

See that $\eqref{eq:avn} - \eqref{eq:lambdavn}$

$$
\boldsymbol{0} = b_1(\lambda_1 - \lambda_n)\boldsymbol{v}_1 + \cdots + b\_{n-1}(\lambda\_{n-1} - \lambda_n)\boldsymbol{v}\_{n-1}
$$

Because the eigenvalues are assumed to be distinct, $\lambda_i - \lambda_n \neq 0$. We also assumed that not all $b_i$'s are zero. This means $\boldsymbol{v}_1, \cdots, \boldsymbol{v}\_{n-1}$ are linearly independent.

If we continue this process, we can eliminate $\boldsymbol{v}\_{n-1}$ from the set, then $\boldsymbol{v}\_{n-2}$, etc. In the end, we'd be able to find that $\boldsymbol{v}_1 = \boldsymbol{0}$, which contradicts that fact that it's an eigenvector. Thus, the $\boldsymbol{v}_i$'s must be linearly independent.

### Theorem 2

{{<hl>}}Eigenvectors of a symmetric matrix (associated with distinct eigenvalues) are orthogonal.{{</hl>}}

**Proof**: for $\boldsymbol{Ax}_1 = \lambda_1\boldsymbol{x}_1$ and $\boldsymbol{Ax}_2 = \lambda_2\boldsymbol{x}_2$, we need to show that if $\lambda_1 \neq \lambda_2$, then $\boldsymbol{x}_1'\boldsymbol{x}_2 = 0$. See that

$$
\begin{gathered}
    \boldsymbol{x}_2'\boldsymbol{Ax}_1 = \lambda_1 \boldsymbol{x}_2'\boldsymbol{x}_1 \\\\
    \boldsymbol{x}_1'\boldsymbol{Ax}_2 = \lambda_2 \boldsymbol{x}_1'\boldsymbol{x}_2
\end{gathered}
$$

The LHS are equal because $\boldsymbol{A}$ is symmetric, so

$$
\lambda_1 \boldsymbol{x}_2'\boldsymbol{x}_1 = \lambda_2 \boldsymbol{x}_1'\boldsymbol{x}_2 = \lambda_2\boldsymbol{x}_2'\boldsymbol{x}_1 \Rightarrow (\lambda_1 - \lambda_2)\boldsymbol{x}_2'\boldsymbol{x}_1 = 0
$$

Since $\lambda_1 \neq \lambda_2$, $\boldsymbol{x}_2'\boldsymbol{x}_1 = 0$.

## Spectral decomposition

The `eigendecomposition` or spectral decomposition is the factorization of a matrix into a canonical form, where the matrix is represented in terms of its eigenvalues and eigenvectors.

Let $\boldsymbol{A}$ be an $n \times n$ matrix with $n$ linearly independent eigenvectors $\boldsymbol{q}_i$, then $\boldsymbol{A}$ can be factorized as

$$
\boldsymbol{A} = \boldsymbol{QDQ}'
$$

where $\boldsymbol{Q}$ is an $n \times n$ matrix whose $i$-th column is the $\boldsymbol{q}_i$, and $\boldsymbol{D}$ is a diagonal matrix whose diagonal elements are the corresponding eigenvalues. We can show that the $\boldsymbol{A}$ matrix is expressed as the sum of $n$ rank 1 matrices:

$$
\boldsymbol{A} = \lambda_1 \boldsymbol{q}_1\boldsymbol{q}_1' + \lambda_2 \boldsymbol{q}_2\boldsymbol{q}_2' + \cdots + \lambda_n \boldsymbol{q}_n\boldsymbol{q}_n' = \sum\_{i=1}^n \lambda_i \boldsymbol{q}_i\boldsymbol{q}_i'
$$

If $\boldsymbol{A}$ has rank $r$, i.e. the first $r$ eigenvalues are non-zero, then $\boldsymbol{A}$ is the sum of the first $r$ rank 1 matrices.

### Low-rank approximation

Suppose $\boldsymbol{A}$ is non-negative definite so that $\lambda_i \geq 0$. The `leading eigenvalues` are the first $R$ eigenvalues that are much larger[^much-larger] than the others:

$$
\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_R \gg \lambda\_{R+1} \geq \cdots \geq \lambda_r > 0
$$

[^much-larger]: The "much larger" here is very subjective, but in real data it happens quite often.

If this is the case, then $\boldsymbol{A}$ can be approximated by the sum of the first $R$ rank 1 matrices $\boldsymbol{A} \approx \sum\_{i=1}^R \lambda_i \boldsymbol{q}_i\boldsymbol{q}_i'$ where $R \ll r$. This is called `low-rank approximation` and is frequently used in machine learning.

### Relationship with previous concepts

{{<alert info>}}
A quick recap on the relationship between eigenvalues and the definiteness of matrices:

-   If all $\lambda_i \geq 0$, then $\boldsymbol{A}$ is non-negative definite.
-   If all $\lambda_i > 0$, then $\boldsymbol{A}$ is positive definite.
-   If all $\lambda_i < 0$, then $\boldsymbol{A}$ is negative definite.
-   If some $\lambda_i > 0$ and some $\lambda_i < 0$, then $\boldsymbol{A}$ is indefinite.

{{</alert>}}

1. Finding the **trace** (again) with the eigendecomposition:

    $$
    tr(\boldsymbol{A}) = tr(\boldsymbol{QDQ}') = tr(\boldsymbol{DQ}'\boldsymbol{Q}) = tr(\boldsymbol{DI}) = tr(\boldsymbol{D}) = \sum\_{i=1}^n \lambda_i
    $$

2. **Determinant**:

    $$
    \begin{aligned}
        \det(\boldsymbol{A}) = \det(\boldsymbol{QDQ}') &= \det(\boldsymbol{Q}) \det(\boldsymbol{D}) \det(\boldsymbol{Q}') \\\\
        &= \det(\boldsymbol{Q}) \det(\boldsymbol{Q}') \det(\boldsymbol{D}) \\\\
        &= \det(\boldsymbol{QQ}') \det(\boldsymbol{D}) \\\\
        &= \det(\boldsymbol{D}) = \prod\_{i=1}^n \lambda_i
    \end{aligned}
    $$

3. Where does the eigendecomposition come from? we can post-multiply $\boldsymbol{Q}$ to both sides:

    $$
    \boldsymbol{AQ} = \boldsymbol{QDQ}'\boldsymbol{Q} = \boldsymbol{QD}
    $$

    Now if we express $\boldsymbol{Q}$ using its column vectors and $\boldsymbol{D}$ its diagonals,

    $$
    \begin{gathered}
        \boldsymbol{A}[\boldsymbol{q}_1, \cdots, \boldsymbol{q}_n] = [\boldsymbol{q}_1, \cdots, \boldsymbol{q}_n] \begin{pmatrix}
            \lambda_1 & & \\\\
            & \ddots & \\\\
            & & \lambda_n
        \end{pmatrix} \\\\
        [\boldsymbol{Aq}_1, \cdots, \boldsymbol{Aq}_n] = [\lambda_1\boldsymbol{q}_1, \cdots, \lambda_n\boldsymbol{q}_n] \\\\
        \boldsymbol{Aq}_i = \lambda_i \boldsymbol{q}_i
    \end{gathered}
    $$

    which is just the eigenvalue-eigenvector pair coming from Eq.$\eqref{eq:ev}$.

4. **Powers**: $\lambda_i^m$ is the eigenvalue for $\boldsymbol{A}^m$.

    $$
    \begin{aligned}
        \boldsymbol{A} &= \boldsymbol{QDQ}' \\\\
        \boldsymbol{A}^2 &= \boldsymbol{QDQ}'\boldsymbol{QDQ}' = \boldsymbol{QD}^2\boldsymbol{Q}' \\\\
        &\quad\vdots \\\\
        \boldsymbol{A}^m &= \boldsymbol{QD}^m\boldsymbol{Q}' \\\\
    \end{aligned}
    $$

    This means for eigenvalues that are smaller than one, if we keep multiplying the matrix by itself then the eigenvalue goes to zero.

5. **Inverse**: if $\boldsymbol{A}$ is non-singular, i.e. none of the eigenvalues are zero,

    $$
    \boldsymbol{A}^{-1} = \left(\boldsymbol{QDQ}'\right)^{-1} = \left(\boldsymbol{Q}'\right)^{-1} \boldsymbol{D}^{-1} \boldsymbol{Q}^{-1} = \boldsymbol{QD}^{-1}\boldsymbol{Q}'
    $$

6. **Square root of a matrix**. The definition is $\boldsymbol{A}^\frac{1}{2}\boldsymbol{A}^\frac{1}{2} = \boldsymbol{A}$. If $\boldsymbol{A}$ is non-negative definite, then

    $$
    \boldsymbol{A}^\frac{1}{2} \triangleq \boldsymbol{QD}^\frac{1}{2}\boldsymbol{Q}', \quad \boldsymbol{D}^\frac{1}{2} = \begin{pmatrix}
        \sqrt{\lambda_1} & & \\\\
        & \ddots & \\\\
        & & \sqrt{\lambda_n}
    \end{pmatrix}
    $$

7. For an **idempotent and symmetric matrix** $\boldsymbol{A}$, if $\lambda$ is an eigenvalue of $\boldsymbol{A}$, $\boldsymbol{Ax} = \lambda\boldsymbol{x}$ for some non-zero vector $\boldsymbol{x}$.

    $$
    \lambda\boldsymbol{x} = \boldsymbol{Ax} = \boldsymbol{A}^2\boldsymbol{x} = \lambda\boldsymbol{Ax} = \lambda^2\boldsymbol{x}
    $$

    So from $\lambda^2 = \lambda$, $\lambda$ can only be 0 or 1. If all $\lambda_i = 1$, $\boldsymbol{A} = \boldsymbol{I}$. Another interesting property is

    $$
    \\#(\lambda_i = 1) = \sum\_{i=1}^n \lambda_i = tr(\boldsymbol{A}) = r(\boldsymbol{A})
    $$

8. For an **orthogonal matrix** $\boldsymbol{P}$, take the norm square:

    $$
    \begin{gathered}
        \boldsymbol{x}'\boldsymbol{P}'\boldsymbol{Px} = \lambda^2\boldsymbol{x}'\boldsymbol{x} \\\\
        \boldsymbol{x}'\boldsymbol{x} = \lambda^2\boldsymbol{x}'\boldsymbol{x} \\\\
        \lambda^2 = 1 \Rightarrow \lambda = \pm 1
    \end{gathered}
    $$

9. We've established that for a symmetric matrix $\boldsymbol{A}$,

    $$
    \begin{aligned}
        \boldsymbol{A} &= \boldsymbol{QDQ}' \\\\
        &= d_1\boldsymbol{q}_1\boldsymbol{q}_1' + \cdots + d_n\boldsymbol{q}_n\boldsymbol{q}_n' \\\\
        &= \sum\_{i=1}^n d_i\boldsymbol{q}_i\boldsymbol{q}_i'
    \end{aligned}
    $$

    Note that each $\boldsymbol{q}_i\boldsymbol{q}_i'$ is an $n \times n$ matrix. If we look at the $j$-th diagonal element of matrix $\boldsymbol{A}$,

    $$
    \boldsymbol{a}\_{jj} = \sum\_{i=1}^n d_i q\_{ji} q\_{ji} = \sum\_{i=1}^n d_i q\_{ji}^2
    $$

    We know that the $q\_{ji}^2$ part is always non-negative, so the diagonal element's sign depends on the corresponding eigenvalue. For off-diagonal elements of $\boldsymbol{A}$,

    $$
    \boldsymbol{a}\_{jk} = TODO
    $$

10. If $\boldsymbol{A}$ has rank $r$ and $r < n$, then we can partition $\boldsymbol{Q}$ into an $n \times r$ matrix $\boldsymbol{Q}_1$ and $n \times (n-r)$ matrix $\boldsymbol{Q}_2$:

    $$
    \begin{aligned}
        \boldsymbol{A} &= \boldsymbol{QDQ}' \\\\
        &= [\boldsymbol{Q}_1, \boldsymbol{Q}_2] \begin{pmatrix}
        \lambda_1 & & & & & \\\\
        & \ddots & & & & \\\\
        & & \lambda_r & & & \\\\
        & & & 0 & & \\\\
        & & & & \ddots & \\\\
        & & & & & 0
    \end{pmatrix} \begin{pmatrix}
        \boldsymbol{Q}_1' \\\\ \boldsymbol{Q}_2'
    \end{pmatrix} \\\\
    &= \boldsymbol{Q}_1 \boldsymbol{D}_1 \boldsymbol{Q}_1' + \boldsymbol{Q}_2 \boldsymbol{0} \boldsymbol{Q}_2'
    \end{aligned}
    $$

    The columns of $\boldsymbol{Q}_1$ are in the column space of $\boldsymbol{A}$, which means $\boldsymbol{Q}_1$ has an orthonormal basis for $\mathcal{C}(\boldsymbol{A})$. Following the same logic, we have $\boldsymbol{Q}_2$ has an orthonormal basis for $\mathcal{N}(\boldsymbol{A})$.

    Thus, {{<hl>}}$\boldsymbol{Q}_1\boldsymbol{Q}_1'$ is the projection matrix to $\mathcal{C}(\boldsymbol{A})$, and $\boldsymbol{Q}_2\boldsymbol{Q}_2'$ is the projection matrix to $\mathcal{N}(\boldsymbol{A})$.{{</hl>}}

11. If $\boldsymbol{A}$ has ordered eigenvalues $\lambda_1 \geq \cdots \geq \lambda_n$, then $\lambda_1 \boldsymbol{I} - \boldsymbol{A}$ is always non-negative definite.

    We can prove this by changing the identity matrix to $\boldsymbol{QQ}'$, where $\boldsymbol{Q}$ is the eigenvector matrix for $\boldsymbol{A}$:

    $$
    \begin{aligned}
        \lambda_1 \boldsymbol{QQ}' - \boldsymbol{QDQ}' &= \boldsymbol{Q}(\lambda_1 \boldsymbol{I}) \boldsymbol{Q}' - \boldsymbol{QDQ}' \\\\
        &= \boldsymbol{Q} \begin{pmatrix}
            \lambda_1 - \lambda_1 & & & \\\\
            & \lambda_1 - \lambda_2 & & \\\\
            & & \ddots & \\\\
            & & & \lambda_1 - \lambda_n
        \end{pmatrix} \boldsymbol{Q}'
    \end{aligned}
    $$

    So $\lambda_1 - \lambda_1, \cdots, \lambda_1 - \lambda_n$ are eigenvalues of $\lambda_1 \boldsymbol{I} - \boldsymbol{A}$. Since all eigenvalues are non-negative, the matrix is non-negative definite.

    Similarly, $\boldsymbol{A} - \lambda_n \boldsymbol{I}$ is non-negative definite.

## Singular value decomposition

The eigendecomposition is a powerful factorization tool, but it only works on square matrices. The <abbr title="Singular Value Decomposition">SVD</abbr> generalizes the eigendecomposition to **any** $m \times n$ matrix.

Suppose $\boldsymbol{A}$ is an $m \times n$ matrix with rank $r$. Obviously we cannot apply eigenvalue decomposition directly on this matrix, but we can try to make it into a square (and symmetric) matrix:

1. $\boldsymbol{AA}'$ is a non-negative definite, symmetric $m \times m$ matrix with rank $r$. Note that $\mathcal{C}(\boldsymbol{AA}') = \mathcal{C}(\boldsymbol{A})$. By eigendecomposition, we have

    $$
    \begin{equation}\label{eq:svd-P-matrix}
        \boldsymbol{AA}' = \boldsymbol{P}\_{m \times m} \boldsymbol{D}\_{m \times m}^* \boldsymbol{P}\_{m \times m}'
    \end{equation}
    $$

    Given it has rank $r$, we can further partition $\boldsymbol{P}$ into $m \times r$ matrix $\boldsymbol{P}_1$ and $m \times (m-r)$ matrix $\boldsymbol{P}_2$:

    $$
    \eqref{eq:svd-P-matrix} = [\boldsymbol{P}_1, \boldsymbol{P}_2] \begin{pmatrix}
        (\boldsymbol{D}_1)\_{r \times r}^* & \boldsymbol{0} \\\\
        \boldsymbol{0} & \boldsymbol{0}
    \end{pmatrix}
    \begin{pmatrix}
        \boldsymbol{P}_1' \\\\ \boldsymbol{P}_2'
    \end{pmatrix} = \boldsymbol{P}_1 \boldsymbol{D}_1^* \boldsymbol{P}_1'
    $$

    The columns of $\boldsymbol{P}_1$ form a basis for $\mathcal{C}(\boldsymbol{A})$.

2. $\boldsymbol{A}'\boldsymbol{A}$ is also non-negative definite, symmetric with rank $r$, but its dimensions are $n \times n$ and we have $\mathcal{R}(\boldsymbol{A}'\boldsymbol{A}) = \mathcal{R}(\boldsymbol{A})$.

    $$
    \begin{equation}\label{eq:svd-Q-matrix}
    \boldsymbol{A}'\boldsymbol{A} = \boldsymbol{Q}\_{n \times n} \boldsymbol{D}^{\*\*} \boldsymbol{Q}' = [\boldsymbol{Q}_1, \boldsymbol{Q}_2] \begin{pmatrix}
            \boldsymbol{D}_1^{\*\*} & \boldsymbol{0} \\\\
            \boldsymbol{0} & \boldsymbol{0}
        \end{pmatrix}
        \begin{pmatrix}
            \boldsymbol{Q}_1' \\\\ \boldsymbol{Q}_2'
        \end{pmatrix} = \boldsymbol{Q}_1 \boldsymbol{D}_1^{\*\*} \boldsymbol{Q}_1'
    \end{equation}
    $$

    Similarly, the $r$ columns of $\boldsymbol{Q}_1$ are a basis for $\mathcal{R}(\boldsymbol{A})$.

Now we want to show that $\boldsymbol{AA}'$ and $\boldsymbol{A}'\boldsymbol{A}$ have the same eigenvalues, i.e. $\boldsymbol{D}_1^* = \boldsymbol{D}_1^{\*\*}$. If $\lambda$ is an eigenvalue of $\boldsymbol{AA}'$,

$$
\boldsymbol{AA}' \boldsymbol{x} = \lambda \boldsymbol{x}
$$

for some $\boldsymbol{x} \neq \boldsymbol{0}$. We're going to pre-multiply $\boldsymbol{A}'$ to both sides:

$$
\begin{gathered}
    \boldsymbol{A}'\boldsymbol{AA}' \boldsymbol{x} = \lambda \boldsymbol{A}'\boldsymbol{x} \\\\
    \boldsymbol{A}'\boldsymbol{Ay} = \lambda\boldsymbol{y}, \quad \boldsymbol{y} = \boldsymbol{A}'\boldsymbol{x}
\end{gathered}
$$

Thus, $\lambda$ is an eigenvalue of $\boldsymbol{A}'\boldsymbol{A}$, QED.

### Definition

The `singular value decomposition` of $\boldsymbol{A}$ decomposes the $m \times n$ matrix into

$$
\boldsymbol{A} = \boldsymbol{PDQ}'
$$

where $\boldsymbol{P}: m \times m$ and $\boldsymbol{Q}: n \times n$ are orthogonal matrices, and $\boldsymbol{D}: m \times n$ is rectangular diagonal[^rectangular-diagonal]. $\boldsymbol{P}$ is called the `left singular vector matrix`, and $\boldsymbol{Q}$ is the `right singular vector matrix`.

[^rectangular-diagonal]: All entries not of the form $d\_{ii}$ are zeros. Think of this as a square diagonal matrix horizontally or vertically concatenated with a zero matrix.

If the rank of $\boldsymbol{A}$ is $r$, then we can partition the orthogonal matrices similarly:

$$
\boldsymbol{A} = [\boldsymbol{P}_1, \boldsymbol{P}_2] \begin{pmatrix}
    \boldsymbol{D}_1 & \boldsymbol{0} \\\\
    \boldsymbol{0} & \boldsymbol{0}
\end{pmatrix} \begin{pmatrix}
    \boldsymbol{Q}_1' \\\\ \boldsymbol{Q}_2'
\end{pmatrix} = \boldsymbol{P}_1 \boldsymbol{D}_1 \boldsymbol{Q}_1'
$$

This is [compact SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition#Reduced_SVDs) and is more economical than the full SVD. The dimensions should be $m \times r$ for $\boldsymbol{P}_1$, $r \times r$ for $\boldsymbol{D}_1$, and $r \times n$ for $\boldsymbol{Q}_1'$. Diagonal elements of $\boldsymbol{D}_1$ are called `singular values`.

Note that $\boldsymbol{P}$ is nothing but the eigenvector matrix of $\boldsymbol{AA}'$, which is exactly the same $\boldsymbol{P}$ matrix in Eq.$\eqref{eq:svd-P-matrix}$. Similarly, $\boldsymbol{Q}$ is the same matrix from Eq.$\eqref{eq:svd-Q-matrix}$, the eigenvector matrix of $\boldsymbol{A}'\boldsymbol{A}$.

We also have $\boldsymbol{D}_1 = \sqrt{\boldsymbol{D}_1^*} = \sqrt{\boldsymbol{D}_1^{\*\*}}$. The singular values of $\boldsymbol{A}$ is the same as the square root of the eigenvalues of $\boldsymbol{AA}'$ or $\boldsymbol{A}'\boldsymbol{A}$. This means no matter what form $\boldsymbol{A}$ takes, its **singular values are strictly positive**. This provides a good measure for the "magnitude" of the matrix.

In practice, singular vectors are often only the part of $\boldsymbol{P}$ or $\boldsymbol{Q}$ that correspond to non-zero singular values, i.e. $\boldsymbol{P}_1$ and $\boldsymbol{Q}_1$. $\boldsymbol{P}_1$ has an orthonormal basis for $\mathcal{C}(\boldsymbol{A})$, and $\boldsymbol{Q}_1$ has an orthonormal basis for $\boldsymbol{R}(\boldsymbol{A})$. Similarly, $\boldsymbol{P}_2$ has an orthonormal basis for $\mathcal{C}(\boldsymbol{A})^\perp$, and $\boldsymbol{Q}_2$ for $\mathcal{N}(\boldsymbol{A})$.

With these, we have

$$
\begin{gathered}
    \boldsymbol{P}_1'\boldsymbol{P}_1 = \boldsymbol{I}_r = \boldsymbol{Q}_1'\boldsymbol{Q}_1 \\\\
    \boldsymbol{P}_1'\boldsymbol{P}_2 = \boldsymbol{0}\_{r \times (m - r)} \\\\
    \boldsymbol{Q}_1'\boldsymbol{Q}_2 = \boldsymbol{0}\_{r \times (n - r)} \\\\
\end{gathered}
$$

Finally, if we represent the singular vectors as individual columns,

$$
\boldsymbol{A} = [\boldsymbol{p}_1, \cdots, \boldsymbol{p}_r] \begin{pmatrix}
    d_1 & & \\\\
    & \ddots & \\\\
    & & d_r
\end{pmatrix} \begin{pmatrix}
    \boldsymbol{q}_1' \\\\ \vdots \\\\ \boldsymbol{q}_r'
\end{pmatrix} = \sum\_{i=1}^r d_i \boldsymbol{p}_i \boldsymbol{q}_i'
$$

which is the sum of $r$ rank 1 matrices. We can approximate this with some $R < r$, namely the `truncated SVD`.

### Moore-Penrose inverse

The Moore-Penrose inverse of a matrix is the most widely known generalization of the inverse matrix. For any $m \times n$ matrix $\boldsymbol{A}$ with rank $r$, if $\boldsymbol{G}$ satisfies the following four criteria:

1. $\boldsymbol{AGA} = \boldsymbol{A}$,
2. $\boldsymbol{GAG} = \boldsymbol{G}$,
3. $\boldsymbol{AG}$ is symmetric,
4. $\boldsymbol{GA}$ is symmetric.

$\boldsymbol{G}$ is the `Moore-Penrose inverse`. It is **uniquely defined**, and is often denoted $\boldsymbol{A}^+$.

A simple example is the rank decomposition $\boldsymbol{A} = \boldsymbol{BT}$ where $\boldsymbol{B}: m \times r$ and $\boldsymbol{T}: r \times n$, then

$$
\boldsymbol{G} = T'(TT')^{-1}(B'B)^{-1}B'
$$

is the Moore-Penrose inverse, and it's invariant to the choice of $\boldsymbol{B}$ and $\boldsymbol{T}$. See that $\boldsymbol{G}$ satisfies all four conditions above.

The Moore-Penrose inverse can be computed using the singular value decomposition. If $\boldsymbol{A} = \boldsymbol{PDQ}'$ is the compact SVD of $\boldsymbol{A}$, then

$$
\boldsymbol{A}^+ = \boldsymbol{QD}^{-1}\boldsymbol{P}'
$$

is the Moore-Penrose inverse of $\boldsymbol{A}$.

So why do we need this special case of generalized inverses? Generalized inverses always exist but are not in general unique. The last two conditions results in the uniqueness of the Moore-Penrose inverse, and it's useful for [computing the least squares solution to a system of linear equations](https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Applications).

## Kronecker product

If $\boldsymbol{A}$ is an $m \times n$ matrix and $\boldsymbol{B}$ is a $p \times q$ matrix, then the `Kronecker product` is the $pm \times qn$ block matrix

$$
\boldsymbol{A} \otimes \boldsymbol{B} = \begin{pmatrix}
    a_{11}\boldsymbol{B} & \cdots & a_{1n}\boldsymbol{B} \\\\
    \vdots & \ddots & \vdots \\\\
    a_{m1}\boldsymbol{B} & \cdots & a_{mn}\boldsymbol{B}
\end{pmatrix}
$$

Note that each of the elements here is a $p \times q$ matrix, so this operation is not commutative. In statistical applications, either or both of $\boldsymbol{A}$ and/or $\boldsymbol{B}$ is the identity matrix, matrix of ones, etc.

For example, in an one-way ANOVA model

$$
y_{ij} = \mu + \alpha_i + \epsilon_{ij}, \quad i = 1, 2, \quad j = 1, 2, 3
$$

In matrix notation, we want to express the model as

$$
\boldsymbol{y} = \boldsymbol{X\beta} + \boldsymbol{\epsilon}
$$

where

$$
\begin{gathered}
    \boldsymbol{y} = (y_{11}, y_{12}, y_{13}, y_{21}, y_{22}, y_{23})' \\\\
    \boldsymbol{\beta} = (\mu, \alpha_1, \alpha_2)'
\end{gathered}
$$

To construct the matrix $\boldsymbol{X}$, we can write it out explicitly

$$
\boldsymbol{X} = \begin{pmatrix}
    1 & 1 & 0 \\\\
    1 & 1 & 0 \\\\
    1 & 1 & 0 \\\\
    1 & 0 & 1 \\\\
    1 & 0 & 1 \\\\
    1 & 0 & 1
\end{pmatrix}
$$

which wasn't too hard because the design is simple, but we can see that things could quickly get out of control for more complicated designs. A neater way to express $\boldsymbol{X}$ is using the Kronecker product:

$$
\boldsymbol{X} = \begin{pmatrix}
    \boldsymbol{1}_6 & \boldsymbol{I}_2 \otimes \boldsymbol{1}_3
\end{pmatrix}
$$
