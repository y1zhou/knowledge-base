---
title: "Generalized Inverse"
date: 2020-10-21T12:45:56-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix"] # list of related tags

slug: "linear-algebra-generalized-inverse"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 80 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Pragser Wildsee, Italy" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "PKjeiYAUm0M" # Unsplash ID of the picture
---

In the previous chapter we discussed the inverse for non-singular matrices.

Let $\boldsymbol{A}$ be an $m \times n$ matrix with rank $r$. $\boldsymbol{G}: n \times m$ is the `generalized inverse` of $\boldsymbol{A}$, denoted by $\boldsymbol{A}^{-}$, that satisfies

$$
\boldsymbol{AGA} = \boldsymbol{A}
$$

For example,

$$
\boldsymbol{A} = \begin{pmatrix}
    1 & 3 & 2 \\\\
    2 & 6 & 4
\end{pmatrix}, \quad
\boldsymbol{G}_1 = \begin{pmatrix}
    1 & 0 \\\\
    0 & 0 \\\\
    0 & 0
\end{pmatrix}, \quad
\boldsymbol{G}_2 = \begin{pmatrix}
    -42 & -1 \\\\
    5 & 3 \\\\
    2 & 2
\end{pmatrix}
$$

Obviously the generalized inverse is not unique. When we talk about generalized inverses, we refer to _all the matrices_ that satisfy $\boldsymbol{AGA} = \boldsymbol{A}$, instead of some specific $\boldsymbol{G}$.

## Finding the generalized inverse

We know that for matrix $\boldsymbol{A}$ with rank $r$, it can be decomposed into $\boldsymbol{B}: m \times r$ and $\boldsymbol{T}: r \times n$, both of rank $r$, such that $\boldsymbol{A} = \boldsymbol{BT}$.

Since $\boldsymbol{B}$ has full column rank, it has a left inverse, i.e. there exists $\boldsymbol{L}: r \times m$ such that $\boldsymbol{LB} = \boldsymbol{I}_r$. Since $\boldsymbol{T}$ has full row rank, there exists $\boldsymbol{R}: n \times r$ such that $\boldsymbol{TR} = \boldsymbol{I}_r$.

Let $\boldsymbol{G} = \boldsymbol{RL}$. Now we may check $\boldsymbol{AGA}$:

$$
\boldsymbol{AGA} = (\boldsymbol{BT})\boldsymbol{RL}(\boldsymbol{BT}) = \boldsymbol{BI}_r\boldsymbol{I}_r\boldsymbol{T} = \boldsymbol{BT} = \boldsymbol{A}
$$

## Properties

-   The first property we should know about is its relation to linear systems.

    If $\boldsymbol{Ax} = \boldsymbol{b}$ is consistent, then $\boldsymbol{x} = \boldsymbol{A}^{-}\boldsymbol{b}$ is a solution. Any solution can be written in this form because $\boldsymbol{b} \in \mathcal{C}(\boldsymbol{A})$, and it can be expressed as $\boldsymbol{b} = \boldsymbol{Av}$ for some $\boldsymbol{v}$. Thus,

    $$
    \boldsymbol{AA}^{-}\boldsymbol{b} = \boldsymbol{AA}^{-}\boldsymbol{Av} = \boldsymbol{Av} = \boldsymbol{b}
    $$

-   Left and right inverses are generalized inverses.

    Let $\boldsymbol{A}: m \times n$ have rank $m$. If there exists $\boldsymbol{R}$ such that $\boldsymbol{AR} = \boldsymbol{I}_m$, then $\boldsymbol{ARA} = \boldsymbol{I}_m\boldsymbol{A} = \boldsymbol{A}$. Similarly, we can show the argument for the full column rank case.

-   Similar to the matrix inverse, we have
    -   $(-\boldsymbol{A})^{-} = -\boldsymbol{A}^{-}$.
    -   $(k\boldsymbol{A})^- = \frac{1}{k}\boldsymbol{A}^-$ where $k \neq 0$.
    -   $(\boldsymbol{A}')^- = (\boldsymbol{A}^-)'$. This can be easily shown by transposing both sides of $\boldsymbol{AGA} = \boldsymbol{A}$.
-   For matrices $\boldsymbol{A}$ and $\boldsymbol{B}$,

    $$
    \mathcal{C}(\boldsymbol{B}) \subset \mathcal{C}(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{B} = \boldsymbol{AA}^-\boldsymbol{B}
    $$

    Showing ($\Leftarrow$) is trivial. For the other direction, we know that $\boldsymbol{B}$ can be written as $\boldsymbol{B} = \boldsymbol{AF}$ for some $\boldsymbol{F}$. Then,

    $$
    \boldsymbol{B} = \boldsymbol{AA}^-\boldsymbol{AF} = \boldsymbol{AA}^-\boldsymbol{B}
    $$

    In terms of vectors, this implicates

    $$
    \boldsymbol{x} \in \mathcal{C}(\boldsymbol{A}) \Longleftrightarrow  \boldsymbol{x} = \boldsymbol{Ay} \text{ for some } \boldsymbol{y} \Longleftrightarrow \boldsymbol{x} = \boldsymbol{AA}^-\boldsymbol{x}
    $$

-   A very important property: $\mathcal{C}(\boldsymbol{AA}^-) = \mathcal{C}(\boldsymbol{A})$.

    Showing ($\subset$) is trivial because $\mathcal{C}(\boldsymbol{AB}) \subset \mathcal{C}(\boldsymbol{A})$. For the other direction, see that

    $$
    \mathcal{C}(\boldsymbol{A}) = \mathcal{C}(\boldsymbol{AA}^-\boldsymbol{A}) \subset \mathcal{C}(\boldsymbol{AA}^-)
    $$

-   $r(\boldsymbol{A}^-) \geq r(\boldsymbol{A})$. We know that the rank cannot be gained by matrix multiplication. Multiplying a full-rank matrix is the only way of _not losing_ rank.

    $$
    r(\boldsymbol{A}^-) \geq r(\boldsymbol{AA}^-) = r(\boldsymbol{A})
    $$

-   If $\mathcal{R}(\boldsymbol{B}) \subset \mathcal{R}(\boldsymbol{A})$, and there exists $\boldsymbol{C}$ such that $\mathcal{C}(\boldsymbol{C}) \subset \mathcal{C}(\boldsymbol{A})$, then $\boldsymbol{BA}^-\boldsymbol{C}$ is `invariant` to the choice of the generalized inverse.

    $$
    \begin{gathered}
        \mathcal{R}(\boldsymbol{B}) \subset \mathcal{R}(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{B} = \boldsymbol{MA} \text{ for some } M \\\\
        \mathcal{C}(\boldsymbol{C}) \subset \mathcal{C}(\boldsymbol{A}) \Longleftrightarrow \boldsymbol{C} = \boldsymbol{AF} \text{ for some } F
    \end{gathered}
    $$

    Then we have

    $$
        \boldsymbol{BA}^- \boldsymbol{C} = \boldsymbol{MAA}^- \boldsymbol{AF} = \boldsymbol{MAF},
    $$

    which is not dependent on $\boldsymbol{A}^-$. This property is used in the regression setting of $\boldsymbol{y} = \boldsymbol{X\beta} + \boldsymbol{\epsilon}$, where

    $$
    \boldsymbol{A} = \boldsymbol{X}'\boldsymbol{X}, \quad \boldsymbol{B} = \boldsymbol{X}, \quad \boldsymbol{C} = \boldsymbol{X}'
    $$

    See that $\mathcal{R}(\boldsymbol{B}) \subset \mathcal{R}(\boldsymbol{A})$ and $\mathcal{C}(\boldsymbol{C}) \subset \mathcal{C}(\boldsymbol{A})$. The predicted values are

    $$
    \hat{\boldsymbol{y}} =\boldsymbol{X}(\boldsymbol{X}'\boldsymbol{X})^- \boldsymbol{X}'\boldsymbol{y} = \boldsymbol{BA}^- \boldsymbol{Cy}
    $$

    This means that although there's infinitely number of choices

    > TODO

-   If $\boldsymbol{A}$ has full row rank $m$, then $\boldsymbol{A}'(\boldsymbol{AA}')^{-1}$ is a right inverse.

    If $r(\boldsymbol{A}) = n$, then $(\boldsymbol{A}'\boldsymbol{A})^{-1}$ is a left inverse.

-   $\boldsymbol{A}'(\boldsymbol{AA}')^-$ is a generalized inverse of $\boldsymbol{A}$, and so is $(\boldsymbol{A}'\boldsymbol{A})^- \boldsymbol{A}$.
