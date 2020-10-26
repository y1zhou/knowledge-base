---
title: "Determinant"
date: 2020-10-26T13:25:26-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: ["Linear Algebra", "Matrix"] # list of related tags

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
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

The `determinant` is a concept about a square matrix $\boldsymbol{A}$, denoted $|\boldsymbol{A}|$ or $\det(\boldsymbol{A})$. It's a scalar value that can be computed from the matrix.

For a $2 \times 2$ matrix, the determinant may be defined as

$$
\boldsymbol{A} = \begin{pmatrix}
    a & b \\\\ c & d
\end{pmatrix}, \quad |\boldsymbol{A}| = ad - bc
$$

## Submatrix

Let $\boldsymbol{A}$ be an $n \times n$ matrix. $\boldsymbol{A}_{-i, -j}$ is an $(n-1) \times (n-1)$ matrix by deleting the $i$-th row and $j$-th column from $\boldsymbol{A}$. For example,

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
