---
title: "Linear Dependence and Independence"
date: 2020-08-31T12:33:24-04:00
summary: "A short piece on linearly dependent and independent sets of vectors." # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Correlation"] # list of related tags

slug: "2-linear-dep-and-indep"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 20 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "-uGzWvOYUTM" # Unsplash ID of the picture
---

The following material on linear dependence and independence, although short, is of fundamental importance in this class.

### Definitions

A set (collection) of column vectors $\boldsymbol{a}_1, \cdots, \boldsymbol{a}_k$ of the same size are `linearly dependent` (LD) if there exist real numbers $x_1, \cdots, x_k$, not all zero, such that

$$
\sum_{i=1}^k x_i \boldsymbol{a}_i = \boldsymbol{0}
$$

If no such scalars exist, the set is called `linearly independent` (LIN). The empty set is considered to be linearly independent.

**Example 1:** $\boldsymbol{a}_1 = (1, 1)^\prime$, $\boldsymbol{a}_2 = (1, 2)^\prime$. We set up a linear combination

$$
\begin{gathered}
    x_1\boldsymbol{a}_1 + x_2\boldsymbol{a}_2 = \boldsymbol{0} \\\\
    x_1 \begin{pmatrix}
        1 \\\\ 1
    \end{pmatrix} + x_2 \begin{pmatrix}
        1 \\\\ 2
    \end{pmatrix} = \begin{pmatrix}
        0 \\\\ 0
    \end{pmatrix} \\\\
    \begin{cases}
        x_1 + x_2 = 0 \\\\
        x_1 + 2x_2 = 0
    \end{cases} \Rightarrow x_1 = 0,\, x_2 = 0
\end{gathered}
$$

Since both coefficients are 0, $\boldsymbol{a}_1$ and $\boldsymbol{a}_2$ are LIN.

**Example 2:** $\boldsymbol{a}_1 = (1, 0, -1)^\prime$, $\boldsymbol{a}_2 = (2, 0, 1)^\prime$, $\boldsymbol{a}_3 = (1, 1, 2)^\prime$.

$$
\begin{gathered}
    x_1\boldsymbol{a}_1 + x_2\boldsymbol{a}_2 + x_3\boldsymbol{a}_3 = \boldsymbol{0} \\\\
    \begin{pmatrix}
        x_1 + 2x_2 + x_3 \\\\
        x_3 \\\\
        -x_1 + x_2 + 2x_3
    \end{pmatrix} = \begin{pmatrix}
        0 \\\\ 0 \\\\ 0
    \end{pmatrix} \\\\
    \Rightarrow x_1 = x_2 = x_3 = 0
\end{gathered}
$$

So again the set of vectors are linearly independent.

**Example 3:** $\boldsymbol{a}_1 = (1, -1, 1)^\prime$, $\boldsymbol{a}_2 = (2, 0, 0)^\prime$, $\boldsymbol{a}_3 = (2, 1, -1)^\prime$.

We can express $\boldsymbol{a}_3$ as $-\boldsymbol{a}_1 + \frac{3}{2}\boldsymbol{a}_2$, so the set of vectors are linearly dependent.

### Important facts

**Lemma 1:** Any set of vectors that include a zero vector is linearly dependent. This can be easily shown by setting all other coefficients to 0 and the coefficient for the zero vector to a non-zero value.

**Lemma 2:** A set of vectors $\\{\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n\\}$ is linearly dependent if and only if there exists $\boldsymbol{a}_i$ such that $\boldsymbol{a}_i$ can be expressed as a linear combination of the rest. This is saying that we have some redundancy in the set of vectors.

To prove this, we first show the "if" ($\Leftarrow$) part. Assume that $\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n$ are LD. By definition, there exists $x_1, \cdots, x_n$ that are not all zeros such that

$$
\sum_{i=1}^n x_i \boldsymbol{a}_i = \boldsymbol{0}
$$

Let $x_k$ be a non-zero number among $x_1, \cdots, x_n$. We then have

$$
x_1\boldsymbol{a}_1 + x_2 \boldsymbol{a}_2 + \cdots + x_k\boldsymbol{a}_k + \cdots + x_n\boldsymbol{a}_n = \boldsymbol{0}
$$

If we divide both sides by $x_k$ and isolate the $k^{th}$ term, we get:

$$
\begin{aligned}
    \boldsymbol{a}_k &= -\left( \frac{x_1}{x_k} \boldsymbol{a}_1 + \cdots + \frac{x\_{k-1}}{x_k} \boldsymbol{a}\_{k-1} + \frac{x\_{k+1}}{x_k}\boldsymbol{a}\_{k+1} + \cdots + \frac{x_n}{x_k}\boldsymbol{a}_n \right) \\\\
    &= y_1\boldsymbol{a}_1 + \cdots + y\_{k-1}\boldsymbol{a}\_{k-1} + y\_{k+1}\boldsymbol{a}\_{k+1} + \cdots + y_n\boldsymbol{a}_n, \quad y_i \triangleq -\frac{x_i}{x_k}
\end{aligned}
$$

Now we show the "only if" part. Suppose we can express some $\boldsymbol{a}_k$ as a linear combination of the rest:

$$
\boldsymbol{a}_k = b_1 \boldsymbol{a}_1 + \cdots + b\_{k-1}\boldsymbol{a}\_{k-1} + b\_{k+1}\boldsymbol{a}\_{k+1} + \cdots + b_n\boldsymbol{a}_n
$$

Moving the RHS to the LHS,

$$
-b_1 \boldsymbol{a}_1 - \cdots - b\_{k-1}\boldsymbol{a}\_{k-1} + \boldsymbol{a}_k - b\_{k+1}\boldsymbol{a}\_{k+1} - \cdots - b_n\boldsymbol{a}_n = \boldsymbol{0}
$$

The coefficient for $\boldsymbol{a}_k$ is not zero, thus the set of vectors are LD.

**Lemma 3:** $\\{\boldsymbol{a}_1, \cdots, \boldsymbol{a}_n\\}$ are LIN if we can't express any vector as a linear combination of the others.

**Lemma 4:** Any set containing a linearly dependent subset is linearly dependent.

**Lemma 5:** Any subset of a linearly independent set is linearly independent.
