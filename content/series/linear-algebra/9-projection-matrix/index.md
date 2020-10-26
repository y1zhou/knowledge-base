---
title: "Projection Matrix"
date: 2020-10-23T13:21:22-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix"] # list of related tags

slug: "linear-algebra-projection-matrix"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 90 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

## Idempotent matrix

Suppose $\boldsymbol{A}$ is a square matrix. We call it an `idempotent matrix` if $\boldsymbol{A}^2 = \boldsymbol{A}$. Some properties of idempotent matrices are:

-   $\boldsymbol{A}'$ is also idempotent.

    $$
    (\boldsymbol{A}')^2 = (\boldsymbol{A}^2)' = \boldsymbol{A}'
    $$

-   $\boldsymbol{A}^2$ is also idempotent.

-   $\boldsymbol{A}^k = \boldsymbol{A}$ where $k = 1, 2, 3, \cdots$, i.e. $\boldsymbol{A}$ is a generalized inverse of $\boldsymbol{A}$. This can be shown by induction.
-   $\boldsymbol{I} - \boldsymbol{A}$ is also idempotent.
-   The only non-singular idempotent matrix is the identity matrix $\boldsymbol{I}$. In general, idempotent matrices are singular.
    -   Some examples of idempotent matrices are $\frac{1}{n}\boldsymbol{J}$ and $(\boldsymbol{I} - \frac{1}{n}\boldsymbol{J})$.
-   For any $m \times n$ matrix $\boldsymbol{A}$, $\boldsymbol{AA}^-$ and $\boldsymbol{A}^-\boldsymbol{A}$ are idempotent.

    $$
    (\boldsymbol{AA}^-)^2 = (\boldsymbol{AA}^-\boldsymbol{A})\boldsymbol{A}^- = \boldsymbol{AA}^-
    $$

-   If $\boldsymbol{A}$ is idempotent, $r(\boldsymbol{A}) = tr(\boldsymbol{A})$. This implies that the trace of an idempotent matrix must be an integer.

    To prove this, recall that $\boldsymbol{B} = \boldsymbol{C} \Longleftrightarrow \boldsymbol{AB} = \boldsymbol{AC}$ if $\boldsymbol{A}$ has full column rank, and $\boldsymbol{B} = \boldsymbol{C} \Longleftrightarrow \boldsymbol{BA} = \boldsymbol{CA}$ if $\boldsymbol{A}$ has full row rank.

    We may combine the two statements above and get $\boldsymbol{B} = \boldsymbol{C} \Longleftrightarrow \boldsymbol{ABD} = \boldsymbol{ACD}$ if $\boldsymbol{A}$ has full column rank and $\boldsymbol{D}$ has full row rank. Now, given $\boldsymbol{A}^2 = \boldsymbol{A}$, suppose $r(\boldsymbol{A}) = r$. We can express $\boldsymbol{A} = \boldsymbol{BT}$ where $\boldsymbol{B}: n \times r$, $\boldsymbol{T}: r \times n$ and $r(\boldsymbol{B}) = r(\boldsymbol{T}) = r$. Thus,

    $$
    \boldsymbol{A}^2 = \boldsymbol{BTBT} = \boldsymbol{BT} = \boldsymbol{BIT}
    $$

    From the equation above we have $\boldsymbol{TB} = \boldsymbol{I}_r$. So,

    $$
    tr(\boldsymbol{A}) = tr(\boldsymbol{BT}) = tr(\boldsymbol{TB}) = tr(\boldsymbol{I}_r) = r = r(\boldsymbol{A})
    $$

## Projection matrix

Let $\boldsymbol{P}$ be an $n \times n$ matrix. For any vector $\boldsymbol{x} \in \mathbb{R}^n$, $\boldsymbol{Px}$ is a projection to a subspace $\mathcal{V} \in \mathbb{R}^n$, i.e. $\boldsymbol{Px} = p(\boldsymbol{y} \mid \mathcal{V})$.

If $\boldsymbol{X}: m \times n$, how do we find the projection matrix to $\mathcal{C}(\boldsymbol{X})$? For any $\boldsymbol{y} \in \mathbb{R}^m$, let $\boldsymbol{z}$ be the projection of $\boldsymbol{y}$ onto $\mathcal{C}(\boldsymbol{X})$:

$$
\boldsymbol{z} = p(\boldsymbol{y} \mid \mathcal{C}(\boldsymbol{X}))
$$

Since $\boldsymbol{z} \in \mathcal{C}(\boldsymbol{X})$, we can write $\boldsymbol{z} = \boldsymbol{Xb}$ for some $\boldsymbol{b} \in \mathbb{R}^n$. Now we have

$$
\boldsymbol{y} - \boldsymbol{z} = \boldsymbol{y} - \boldsymbol{Xb} \perp \mathcal{C}(\boldsymbol{X})
$$

which means it's orthogonal to every column of $\boldsymbol{X}$, thus

$$
\begin{gathered}
    \boldsymbol{X}'(\boldsymbol{y} - \boldsymbol{Xb}) = \boldsymbol{0} \\\\
    \boldsymbol{X}'\boldsymbol{y} - \boldsymbol{X}'\boldsymbol{Xb} = \boldsymbol{0} \\\\
    \boldsymbol{X}'\boldsymbol{Xb} = \boldsymbol{X}'\boldsymbol{y} \\\\
    \boldsymbol{b} = (\boldsymbol{X}'\boldsymbol{X})^-\boldsymbol{X}'\boldsymbol{y}
\end{gathered}
$$

So $\boldsymbol{z} = \boldsymbol{Xb} = \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^-\boldsymbol{X}'\boldsymbol{y}$ is the projection of $\boldsymbol{y}$ onto $\mathcal{C}(\boldsymbol{X})$, and $\boldsymbol{z} = \boldsymbol{Py}$ where

$$
P = \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^-\boldsymbol{X}'
$$

is the projection matrix to $\mathcal{C}(\boldsymbol{X})$.

If $r(\boldsymbol{X}'\boldsymbol{X}) = n$, $\boldsymbol{X}$ has full column rank $n$ and we can use the matrix inverse for the projection matrix

$$
P = \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^{-1}\boldsymbol{X}'
$$

No matter which inverse we end up using, $\boldsymbol{P}$ is unique and is invariant of the choice of the generalized inverse. We can prove this using the 7th fact in [the generalized inverse chapter]({{< relref "../8-generalized-inverse/index.md#properties" >}}), which can help us show that

$$
\mathcal{R}(\boldsymbol{X}) \subset \mathcal{R}(\boldsymbol{X}'\boldsymbol{X}), \quad \mathcal{C}(\boldsymbol{X}') \subset \mathcal{C}(\boldsymbol{X}'\boldsymbol{X})
$$

### Properties

Let $\boldsymbol{P}: n \times n$ be the projection matrix to $\mathcal{C}(\boldsymbol{X})$ where $\boldsymbol{X}: n \times p$.

1. $\boldsymbol{PX} = \boldsymbol{X}$.

    $$
    \begin{aligned}
        \boldsymbol{PX} &= \boldsymbol{P} [\boldsymbol{x}_1, \cdots, \boldsymbol{x}_p] \\\\
        &= [\boldsymbol{Px}_1, \cdots, \boldsymbol{Px}_p] \\\\
        &= [\boldsymbol{x}_1, \cdots, \boldsymbol{x}_p] = \boldsymbol{X}
    \end{aligned}
    $$

    From this fact, we can see that $(\boldsymbol{X}'\boldsymbol{X})^- \boldsymbol{X}'$ is a generalized inverse of $\boldsymbol{X}$.

2. $\boldsymbol{P}$ is symmetric. Note that $(\boldsymbol{X}'\boldsymbol{X})^-$ may not be symmetric.

    $$
    \begin{aligned}
        \boldsymbol{P}' &= \left( \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^-\boldsymbol{X}' \right)' \\\\
        &= \boldsymbol{X}\left((\boldsymbol{X}'\boldsymbol{X})^-\right)'\boldsymbol{X}'
    \end{aligned}
    $$

    Since $\left((\boldsymbol{X}'\boldsymbol{X})^-\right)'$ is a generalized inverse of $(\boldsymbol{X}'\boldsymbol{X})' = \boldsymbol{X}'\boldsymbol{X}$, $\boldsymbol{P}'$ is a projection matrix to $\mathcal{C}(\boldsymbol{X})$. As the projection matrix is unique, $\boldsymbol{P} = \boldsymbol{P}'$.

3. $\boldsymbol{X}'\boldsymbol{P} = \boldsymbol{X}'$. Take the transpose of (1) and use (2).
4. $\boldsymbol{P}$ is idempotent. We can use (1) to show this:

    $$
    \begin{aligned}
        \boldsymbol{P}^2 &= \boldsymbol{P} \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^- \boldsymbol{X}' \\\\
        &= \boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^- \boldsymbol{X}' = \boldsymbol{P}
    \end{aligned}
    $$

5. $\mathcal{C}(\boldsymbol{P}) = \mathcal{C}(\boldsymbol{X})$. By the form of $\boldsymbol{P}$ we can obviously see that $\mathcal{C}(\boldsymbol{P}) \subset \mathcal{C}(\boldsymbol{X})$. $\mathcal{C}(\boldsymbol{X}) \subset \mathcal{C}(\boldsymbol{P})$ can be seen from $\mathcal{C}(\boldsymbol{X}) = \mathcal{C}(\boldsymbol{PX}) \subset \mathcal{C}(\boldsymbol{P})$.
6. $r(\boldsymbol{P}) = r(\boldsymbol{X})$. Since $\boldsymbol{P}$ is idempotent, $r(\boldsymbol{P}) = tr(\boldsymbol{P}) = r(\boldsymbol{X})$.
7. $\boldsymbol{I} - \boldsymbol{P}$ is the projection matrix to $\mathcal{C}(\boldsymbol{X})^\perp$.

    For any $\boldsymbol{y} \in \mathbb{R}^n$, we have

    $$
    \begin{aligned}
        \boldsymbol{y} &= p(\boldsymbol{y} \mid \mathcal{C}(\boldsymbol{X})) + p(\boldsymbol{y} \mid \mathcal{C}(\boldsymbol{X})^\perp) \\\\
        &= \boldsymbol{Py} + (\boldsymbol{I} - \boldsymbol{P})\boldsymbol{y}
    \end{aligned}
    $$

8. If $\boldsymbol{A}$ is symmetric and idempotent, then $\boldsymbol{A}$ is a projection matrix to $\mathcal{C}(\boldsymbol{A})$ and vice versa.

    We've already shown ($\Leftarrow$) before. To show ($\Rightarrow$), see that the projection matrix to $\mathcal{C}(\boldsymbol{A})$ is

    $$
    \begin{aligned}
        \boldsymbol{A}(\boldsymbol{A}'\boldsymbol{A})^-\boldsymbol{A}' &= \boldsymbol{A}(\boldsymbol{A}^2)^- \boldsymbol{A} \\\\
        &= \boldsymbol{AA}^- \boldsymbol{A} \\\\
        &= \boldsymbol{A}
    \end{aligned}
    $$

9. If $\boldsymbol{X}$ is non-singular (or of full row rank $n$), then $\boldsymbol{P} = \boldsymbol{I}_n$. This is because $\mathcal{C}(\boldsymbol{X}) = \mathbb{R}^n$ and the projection stays the same.
