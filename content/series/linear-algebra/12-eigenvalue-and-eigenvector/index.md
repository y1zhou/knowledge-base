---
title: "Eigenvalues and Eigenvectors"
date: 2020-11-18T18:45:52-05:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "linear-algebra-eigenvalues-and-eigenvectors"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 120 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
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

So $\boldsymbol{x} = (1, 0)'$ is our eigenvector, where $x_1 = 1$ because we want to make the norm 1. Sometimes we would find more than one (say, $k$) eigenvector associated with the same eigenvalue, then we would say the eigenvalue has geometric multiplicity of $k$. The `geometric multiplicity` is the dimension of the nullspace of $\boldsymbol{A} - \lambda\boldsymbol{I}$.

In general, the algebraic multiplicity and geometric multiplicity of an eigenvalue can differ, but the geometric multiplicity can **never exceed** the algebraic multiplicity.
