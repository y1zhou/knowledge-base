---
title: "Orthogonalization"
date: 2020-09-28T21:55:43-04:00
summary: "Introducing the Gram-Schmidt process, a method for constructing an orthogonal basis given a non-orthogonal basis." # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Operations"] # list of related tags

slug: "linear-algebra-orthogonalization"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 43 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Corner of a building." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "ggKBImNKmiU" # Unsplash ID of the picture
---

## Gram-Schmidt orthogonalization

The `Gram-Schmidt process` is a method for orthogonalizing a set of vectors. We start with a non-orthogonal basis $\\{\boldsymbol{u}_1, \cdots, \boldsymbol{u}_k\\}$ for $V$, and our goal is to find $\\{\boldsymbol{v}_1, \cdots, \boldsymbol{v}_k\\}$, an orthogonal basis for $V$.

1. The first member in the orthogonal basis is the same as the first member[^orthogonal-basis] in the given basis: $\boldsymbol{v}_1 \triangleq \boldsymbol{u}_1$.
2. We need $\boldsymbol{v}_2$ such that $\boldsymbol{v}_2 \perp \boldsymbol{v}_1$, and $\mathcal{L}(\boldsymbol{u}_1, \boldsymbol{u}_2) = \mathcal{L}(\boldsymbol{v}_1, \boldsymbol{v}_2)$. We can show that $\boldsymbol{v}_2$ can be found by taking the difference between $\boldsymbol{u}_2$ and its projection on $\boldsymbol{v}_1$:

    $$
    \begin{aligned}
        \boldsymbol{v}_2 &\triangleq \boldsymbol{u}_2 - P(\boldsymbol{u}_2 \mid \boldsymbol{v}_1) \\\\
        &= \boldsymbol{u}_2 - \left( \frac{\boldsymbol{u}_2 \cdot \boldsymbol{v}_1}{\\|\boldsymbol{v}_1\\|^2} \right)\boldsymbol{v}_1
    \end{aligned}
    $$

3. We find $\boldsymbol{v}_3$ such that $\boldsymbol{v}_3 \perp \mathcal{L}(\boldsymbol{v}_1, \boldsymbol{v}_2)$ and $\mathcal{L}(\boldsymbol{u}_1, \boldsymbol{u}_2, \boldsymbol{u}_3) = \mathcal{L}(\boldsymbol{v}_1, \boldsymbol{v}_2, \boldsymbol{v}_3)$. Similar to the last step,

    $$
    \begin{aligned}
        \boldsymbol{v}_3 &\triangleq \boldsymbol{u}_3 - P(\boldsymbol{u}_3 \mid \mathcal{L}(\boldsymbol{v}_1, \boldsymbol{v}_2)) \\\\
        &= \boldsymbol{u}_3 - \left( \frac{\boldsymbol{u}_3 \cdot \boldsymbol{v}_1}{\\|\boldsymbol{v}_1\\|^2} \right)\boldsymbol{v}_1 - \left( \frac{\boldsymbol{u}_3 \cdot \boldsymbol{v}_2}{\\|\boldsymbol{v}_2\\|^2} \right)\boldsymbol{v}_2
    \end{aligned}
    $$

4. Continue this until we find

    $$
    \boldsymbol{v}\_k = \boldsymbol{u}\_k - P(\boldsymbol{u}\_k \mid \mathcal{L}(\boldsymbol{v}_1, \cdots, \boldsymbol{v}\_{k-1}))
    $$

[^orthogonal-basis]: It doesn't matter which one you start with, but the result may look different when different vectors are chosen.

For example, say we have a linearly independent set of vectors

$$
\begin{gathered}
    \boldsymbol{u}_1 = (1, 1, 1, 1)^\prime \\\\
    \boldsymbol{u}_2 = (1, 1, 5, 5)^\prime \\\\
    \boldsymbol{u}_3 = (4, 0, 12, 8)^\prime
\end{gathered}
$$

and we want to find an orthogonal basis $V = \mathcal{L}(\boldsymbol{u}_1, \boldsymbol{u}_2, \boldsymbol{u}_3) \subset \mathbb{R}^4$.

$$
\begin{aligned}
    \boldsymbol{v}_1 &= (1, 1, 1, 1)^\prime \\\\
    \boldsymbol{v}_2 &= (1, 1, 5, 5)^\prime - \frac{12}{4}(1, 1, 1, 1)^\prime = (-2, -2, 2, 2)^\prime \\\\
    \boldsymbol{v}_3 &= (4, 0, 12, 8)^\prime - \frac{4 + 12 + 8}{4}(1, 1, 1, 1)^\prime - \frac{-8 + 24 + 16}{16}(-2, -2, 2, 2)^\prime = \cdots
\end{aligned}
$$

## Orthogonal complement of subspace

Suppose $V$ is a subspace of $\mathbb{R}^n$. The `orthogonal complement` of $V$ is defined as a set of all vectors that are orthogonal to $V$, often read "$V$ perp":

$$
V^\perp \triangleq \\{\boldsymbol{x}: \boldsymbol{x} \in \mathbb{R}^n, \boldsymbol{x} \perp V \\}
$$

For example, if $V = \mathcal{L}\left(((1, 0)^\prime \right)$, $V^\perp = \mathcal{L}\left((0, 1)^\prime \right)$.

$V \cup V^\perp$ is **not** a subspace, but is a set in $\mathbb{R}^2$. $(3, 4)^\prime \in \mathbb{R}^2$ does not belong to $V \cup V^\perp$, but we can express it as

$$
(3, 4)^\prime = \underbrace{(3, 0)^\prime}\_{\in V} + \underbrace{(0, 4)^\prime}\_{\in V^\perp}
$$

**Theorem:** $V^\perp$ is a subspace of $\mathbb{R}^n$. To prove this we just need to show that it's closed under linear combination. For any $\boldsymbol{x}_1, \boldsymbol{x}_2 \in V^\perp$,

$$
\alpha_1 \boldsymbol{x}_1 + \alpha_2 \boldsymbol{x}_2 \in V^\perp \quad \forall \alpha_1, \alpha_2 \in \mathbb{R}
$$

because $\alpha_1 \boldsymbol{x}_1 + \alpha_2 \boldsymbol{x}_2 \perp V$. For any $\boldsymbol{y} \in V$,

$$
(\alpha_1 \boldsymbol{x}_1 + \alpha_2 \boldsymbol{x}_2) \cdot \boldsymbol{y} = \alpha_1 \boldsymbol{x}_1 \boldsymbol{y} + \alpha_2 \boldsymbol{x}_2 \boldsymbol{y} = 0
$$

**Theorem:** For any vector $\boldsymbol{x} \in \mathbb{R}^n$, there exists $\boldsymbol{u}$ and $\boldsymbol{v}$ such that $\boldsymbol{x} = \boldsymbol{u} + \boldsymbol{v}$ where $\boldsymbol{u} \in V$ and $\boldsymbol{v} \in V^\perp$.

Let $\\{\boldsymbol{u}_1, \cdots, \boldsymbol{u}_k \\}$ be an orthogonal basis of $V$. By Gram-Schmidt, we can find

$$
\\{ \boldsymbol{u}_1, \cdots, \boldsymbol{u}_k, \boldsymbol{u}\_{k+1}, \cdots, \boldsymbol{u}_n \\}
$$

to be an orthogonal basis of $\mathbb{R}^n$. Then, for any $\boldsymbol{x} \in \mathbb{R}^n$, we can write it as a linear combination of the basis vectors:

$$
\boldsymbol{x} = \underbrace{\alpha_1 \boldsymbol{u}_1 + \cdots + \alpha_k \boldsymbol{u}_k}\_{\text{lin. comb. of }\boldsymbol{u}_1, \cdots, \boldsymbol{u}_k \in V} + \underbrace{\alpha\_{k+1} \boldsymbol{u}\_{k+1} + \cdots + \alpha_n \boldsymbol{u}_n}\_{\text{lin. comb. of } \boldsymbol{u}\_{k+1}, \cdots, \boldsymbol{u}_n \in V^\perp}
$$

This implies that $dim(V^\perp) = n-k = n - dim(V)$. We can also see that $(V^\perp)^\perp = V$.

For example, $V = \mathcal{L} \left( \boldsymbol{v}_1 = (1, 1, 1, 1)^\prime, \boldsymbol{v}_2 = (1, 1, -1, -1)^\prime \right)$. We first take two linearly independent vectors that are not in $V$, e.g. $(1, 0, 0, 0)^\prime$ and $(0, 0, 1, 0)^\prime$[^finding-lin-vectors]. Then,

[^finding-lin-vectors]: Tips for constructing this: use a lot of zeros. If you got this step wrong, the last step in G-S would give a vector that's already seen.

$$
\left\\{ (1, 1, 1, 1)^\prime, (1, 1, -1, -1)^\prime, (1, 0, 0, 0)^\prime, (0, 0, 1, 0)^\prime \right\\}
$$

is a basis of $\mathbb{R}^4$. Using G-S to orthogonalize them, we have

$$
\begin{gathered}
    \boldsymbol{v}_3 = \left( 1, -1, 0, 0 \right)^\prime \\\\
    \boldsymbol{v}_4 = \left( 0, 0, 1, -1 \right)^\prime
\end{gathered}
$$

and $V^\perp = \mathcal{L}(\boldsymbol{v}_3, \boldsymbol{v}_4)$. Take

$$
\begin{aligned}
    \boldsymbol{x} &= (1, -2, 3, -4)^\prime \\\\
    &= \boldsymbol{u} + \boldsymbol{v} \quad \boldsymbol{u} \in V, \boldsymbol{v} \in V^\perp \\\\
    &= P(\boldsymbol{x} \mid V) + P(\boldsymbol{x} \mid V^\perp) \\\\
    &= P(\boldsymbol{x} \mid \boldsymbol{v}_1) + P(\boldsymbol{x} \mid \boldsymbol{v}_2) + P(\boldsymbol{x} \mid \boldsymbol{v}_3) + P(\boldsymbol{x} \mid \boldsymbol{v}_4) \\\\
    &= \left( \frac{1-2+3-4}{4}\boldsymbol{v}_1 + \frac{1-2-3+4}{4}\boldsymbol{v}_2 \right) + \left(  \frac{1+2}{2}\boldsymbol{v}_3 + \frac{3 + 4}{2} \boldsymbol{v}_4 \right) \\\\
    &= \left( -\frac{1}{2}(1, 1, 1, 1)^\prime + 0(1, 1, -1, -1)^\prime \right) + \left( \frac{3}{2}(1, -1, 0, 0)^\prime + \frac{7}{2} (0, 0, 1, -1)^\prime \right) \\\\
    &= \left( -\frac{1}{2}, -\frac{1}{2}, -\frac{1}{2}, -\frac{1}{2} \right)^\prime + \left( \frac{3}{2}, -\frac{3}{2}, \frac{7}{2}, -\frac{7}{2} \right)^\prime
\end{aligned}
$$

Note that if we chose different LIN vectors from $(1, 0, 0, 0)^\prime$ and $(0, 0, 1, 0)^\prime$, our $\boldsymbol{v}_3$ and $\boldsymbol{v}_4$ may change, but the two projections are not going to change.
