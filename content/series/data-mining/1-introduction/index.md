---
title: "Overview of the Course"
date: 2021-01-14T15:55:00-04:00
summary: "" # appears in list of posts
categories: ["Data Mining"] # main category; shown in post metadata
tags: ["Data Mining"] # list of related tags

slug: "data-mining-introduction"
toc: false # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 10 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "" # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "" # Unsplash ID of the picture
---

Data mining is largely based on machine learning. A machine learning algorithm learns _patterns_ from a set of provided examples (data). Learning is often spelled out in terms of approximating an ideal target function $t$ with another function $h$. In ML $h$ is often called a `model` or `hypothesis`.

We use machine learning to develop models, and apply what is learned to future examples. Learning is often an iterative process. We begin with some initial hypothesis $h$, and with each iteration, $h$ is modified to become a better approximation of $t$. Learning continues until some stopping criterion is met.

In a typical scenario, we have a set of examples as the training set. Each example has some attributes, and among the attributes is one we are interested in (the target). The other attributes are input features to the function to be learned.

The ML algorithm uses the examples in the training set to learn a function from the inputs to the target. A good function yields the correct value of the target attribute _not just on the training data but on future data sets_. In other words, a good model **generalizes** well.

The function learned could be one of many things -- linear or polynomial equations, decision trees, neural networks, set of rules, etc. Similarly, the input to the function could take different forms. Mostly we will work with tabular data where each row corresponds a combination of input and output, and each column represents an attribute. Often attributes are heterogeneous, meaning that they have different data types and scales.
