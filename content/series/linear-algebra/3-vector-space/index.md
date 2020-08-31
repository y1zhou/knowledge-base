---
title: "Vector Space"
date: 2020-08-31T13:13:34-04:00
summary: "" # appears in list of posts
categories: ["Linear Algebra"] # main category; shown in post metadata
tags: [] # list of related tags

slug: "linear-algebra-vector-space"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: true

weight: 30 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Ocean clouds seen from space." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "yZygONrUBe8" # Unsplash ID of the picture
---

## Definitions

### Vector space

A non-empty set of vectors $V$ is called a `vector space` (or linear space) if:

1. For all $\boldsymbol{x}, \boldsymbol{y} \in V$, $\boldsymbol{x} + \boldsymbol{y} \in V$, i.e. $V$ is closed under addition.
2. For all $\boldsymbol{x} \in V$, $k \in \mathbb{R}$, $k\boldsymbol{x} \in V$, i.e. $V$ is closed under scalar multiplication.

The second condition implies that the zero vector is always included in any vector space $V$. Alternatively, we may say that $V$ is closed under linear combination, i.e.

$$
k_1 \boldsymbol{x}_1 + \cdots + k_n \boldsymbol{x}_n \in V
$$

for all $\boldsymbol{x}_1, \cdots, \boldsymbol{x}_n \in V$ and $k_1, \cdots, k_n \in \mathbb{R}$.

For example, below are some valid/invalid vector spaces:

1. $V = \\{(x, y)^\prime: x \in \mathbb{R}, y \in \mathbb{R}\\}$ is a vector space.
2. $V = \\{(x, y)^\prime: x > y\\}$ is not a vector space.
3. $V = \\{(x, y, 0)^\prime: x \in \mathbb{R}, y \in \mathbb{R}\\}$ is a vector space.
4. $V = \\{(x, y, z)^\prime: 2x = y, x \in \mathbb{R}, y \in \mathbb{R}, z \in \mathbb{R}\\}$

### Subspace

A set $W$ is a `subspace` of a vector space $V$ if $W \subset V$, and $W$ itself is a vector space.

For example, $V = \mathbb{R}^2 = \\{(x, y)^\prime \mid x \in \mathbb{R}, y \in \mathbb{R}\\}$ is a vector space. Let

$$
W_1 = \\{ (x, y)^\prime: y-2x = 0, x, y \in \mathbb{R} \\}
$$
