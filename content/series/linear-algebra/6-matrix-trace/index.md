---
title: "Matrix Trace"
date: 2020-10-07T13:07:16-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "linear-algebra-matrix-trace"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 60 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

The `trace` of a square matrix $\boldsymbol{A}: n \times n$ is the sum of the diagonals:

$$
tr(\boldsymbol{A}) = \sum_{i=1}^n a_{ii}
$$

The trace is related to many important features of matrices, and we're going to discuss some basic properties first.

## Properties

-   **Trace is commutative**, i.e. $tr(\boldsymbol{AB}) = tr(\boldsymbol{BA})$. Note that if $\boldsymbol{A}: m \times n$ and $\boldsymbol{B}: n \times m$, $\boldsymbol{AB}$ would be $m \times m$ and $\boldsymbol{BA}$ would be $n \times n$.

$$
\begin{aligned}
    LHS &= \sum_{i=1}^m (\boldsymbol{AB})_{ii} \\\\
    &= \sum_{i=1}^m \left( \sum_{i=1}^n a_{ij}b_{ji} \right) \\\\
    &= \sum_{j=1}^n \sum_{i=1}^m b_{ji}a_{ij} \\\\
    &= \sum_{j=1}^n (\boldsymbol{BA})_{jj} = tr(\boldsymbol{BA})
\end{aligned}
$$

For more than two matrices, note that $tr(\boldsymbol{ABC}) \neq tr(\boldsymbol{CBA})$. Instead,

$$
tr(\boldsymbol{ABC}) = tr(\boldsymbol{A}(\boldsymbol{BC})) = tr(\boldsymbol{BCA})
$$

Similarly, we have

$$
tr(\boldsymbol{ABCDE}) = tr(\boldsymbol{CDEAB}) = tr(\boldsymbol{EABCD}) = \cdots
$$

-   $tr(\boldsymbol{A}^\prime \boldsymbol{A}) = tr(\boldsymbol{AA}^\prime) \geq 0$.

This is equal to the sum of all entries squared.

-   $tr(\boldsymbol{A}^\prime\boldsymbol{A}) = 0$ is equivalent to $\boldsymbol{A} = 0$.

Showing $\Leftarrow$ is trivial. For $\Rightarrow$, the trace being zero means that the sum of all entries squared of $\boldsymbol{A}$ is zero, thus $a_{ij} = 0$ for all $i, j$.

We can also use this fact to show that $tr(\boldsymbol{A}^\prime \boldsymbol{A}) = 0 \Leftrightarrow \boldsymbol{A}^\prime\boldsymbol{A} = 0$.

-   $\boldsymbol{A}^\prime \boldsymbol{AB} = \boldsymbol{A}^\prime \boldsymbol{AC} \Leftrightarrow \boldsymbol{AB} = \boldsymbol{AC}$.

This one is related to future statistical analyses. Again $\Leftarrow$ is obvious. To prove $\Rightarrow$, we will show $\boldsymbol{AB} = \boldsymbol{AC}$ by seeing that $tr((\boldsymbol{AB} - \boldsymbol{AC})^\prime (\boldsymbol{AB} - \boldsymbol{AC})) = 0$.

We start with $\boldsymbol{A}^\prime \boldsymbol{AB} - \boldsymbol{A}^\prime \boldsymbol{AC} = 0$.

$$
\begin{gathered}
    \boldsymbol{A}^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    (\boldsymbol{B}^\prime - \boldsymbol{C}^\prime) \boldsymbol{A}^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    (\boldsymbol{AB} - \boldsymbol{AC})^\prime (\boldsymbol{AB} - \boldsymbol{AC}) = 0 \\\\
    \boldsymbol{AB} - \boldsymbol{AC} = 0
\end{gathered}
$$

We will use this fact to show that $\mathcal{C}(\boldsymbol{A}) = \mathcal{C}(\boldsymbol{AA}^\prime)$.
