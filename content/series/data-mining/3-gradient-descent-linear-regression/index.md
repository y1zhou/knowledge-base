---
title: "Gradient Descent and Linear Regression"
date: 2021-02-11T15:55:00-04:00
summary: "We implement linear regression using gradient descent, a general optimization technique which in this case can find the global minimum." # appears in list of posts
categories: ["Data Mining"] # main category; shown in post metadata
tags: ["Data Mining", "Regression"] # list of related tags

slug: "data-mining-gradient-descent-linear-regression"
toc: false # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 30 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Grand Canyon National Park, United States." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "1eAbNB7NWEw" # Unsplash ID of the picture
---

In a classification problem, we aim to learn a function that correctly places a newly seen instance into a class based upon the instance's attributes. This is in contrast to `regression`, where we attempt to learn a real-valued function, i.e. `model`, from attributes $x_1, \cdots, x_k$ to a target attribute $y$.

## Linear models

In a linear model, the `hypothesis function` has the form

<div>
$$
    h_\theta (\boldsymbol{x}) = \theta_0 + \sum_{j=1}^k \theta_j x_j,
$$
</div>

where each $\theta_j$ is a weight, and $h_\theta(x)$ (the hypothesis) is the estimated value of $y$. If we define $x_0$ to be a constant 1, then the function can be written more succinctly:

<div>
$$
    h_\theta (\boldsymbol{x}) = \sum_{k=0}^k \theta_j x_j
$$
</div>

The function $h_\theta(x)$ defines a $k$-dimensional `hyperplane`. If $k=1$, it's a line separating a 2d space into two pieces. If $k=2$, then it's a 2d plane, etc. Intuitively, each weight $\theta_j$ represents the "slope" of the plane in the $j$-th dimension.

The problem in constructing a linear model is to find weights that fit the data. Specifically, we want to minimize the squared error between each `$y^{(i)}$` and `$h_\theta (x^{(i)})$`, the actual and predicted values for observation $i$. The term `$y^{(i)} - h_\theta (x^{(i)})$` is called the `residual` of `$y^{(i)}$` with respect to `$h_\theta (x^{(i)})$`. The problem can be cast as minimizing the `cost function` $J(\theta)$:

<div>
$$
    J(\theta) = \frac{1}{2}\sum_{i=1}^n \left( h_\theta (x^{(i)}) - y^{(i)} \right)^2 = \frac{1}{2}\sum_{i=1}^n \left( \sum_{j=0}^k \theta_j x_j^{(i)} - y^{(i)} \right)^2
$$
</div>

The $1/2$ term doesn't change what the best weights are; it ultimately simplifies calculations. Note that $J$ is a function of $\theta$. The space of possible values for $\theta$ defines an `error surface`, and we want to find the values for $\theta$ that correspond to the lowest point on that surface.

## Gradient descent

Minimizing error is a standard objective in machine learning. For linear regression, there exists a unique `global minimum` to $J(\theta)$, and analytic techniques can be used to find it[^linear-regression]. For other functions, techniques such as `gradient descent` are more appropriate.

[^linear-regression]:
    In matrix notation, if our model is $\boldsymbol{y} = \boldsymbol{X\beta} + \boldsymbol{\epsilon}$, then the weight coefficients are

    <div>
    $$
        \boldsymbol{\beta} = (\boldsymbol{X}'\boldsymbol{X})^{-1}\boldsymbol{X}'\boldsymbol{y}
    $$
    </div>

    Note that the first column of $\boldsymbol{X}$ is a column of ones, and $\beta_0$ is the intercept term.

In gradient descent, the gradient vector at a given point on the error surface represents the direction of greatest rate of increase in $J$ at the point:

<div>
$$
    \nabla J = \langle \frac{\partial J}{\partial \theta_0},\ldots, \frac{\partial J}{\partial \theta_k}\rangle
$$
</div>

Here the point is defined by specific values for $\theta_0, \cdots, \theta_k$. Intuitively, `$\frac{\partial J}{\partial \theta_j}$` represents the slope of the surface in the $j$-th dimension.

The idea behind gradient descent is to **adjust the current weights in the opposite direction of the gradient, thereby decreasing the error value**. The weights are initialized (typically to 0), and then iteratively adjusted. The equation to update each existing weight $\theta_j$ to a new value $\theta_j'$ is:

<div>
$$
    \begin{equation}\label{eq:gradient-descent}
        \theta_j' = \theta_j - \alpha \frac{\partial J}{\partial \theta_j} = \theta_j - \alpha \left( \sum_{i=1}^n \left[\left( \sum_{j=0}^k \theta_j x_j^{(i)} - y^{(i)} \right) x_j^{(i)} \right] \right)
    \end{equation}
$$
</div>

where $\alpha$ is a small real-valued constant called the `learning rate`[^learning-rate]. Note that all weights are updated in parallel. If $\nabla J = 0$ at a given point, then further updates yield no improvement. The point constitutes a `local minimum` for $J(\theta)$. Gradient descent stops when the local minimum is reached. In the context of linear regression, the cost function $J$ is `convex`, so there is only one minimum -- the global minimum. In other cases, gradient descent is not guaranteed to find the global minimum.

[^learning-rate]: If $\alpha$ is too large, then $J$ might not converge; it either increases without a bound or oscillates between points. If $\alpha$ is too small, then gradient descent may take forever to converge.

### Batch and stochastic gradient descent

In the simplest form of gradient descent, all instances in teh data are examined before updates are made. This is called `batch gradient descent`. We take the average of the gradients of all instances and use that mean gradient to update our parameters.

A variant is called `stochastic gradient descent` (SGD), where a randomly chosen instance[^mini-batch] is used instead of the entire data set. The parameters are updated using the computed gradients. This often reduces error more quickly. The algorithm might not converge to a minimum, but the results are typically good approximations.

[^mini-batch]: In `mini-batch gradient descent`, the two are combined where we choose a random sample of the training data, and use the mean gradient to update the parameters.

### Scaling

When fitting a function of several attributes, it's often the case that attributes will have very different magnitudes. This can drastically reduce the efficiency of gradient descent. As such, it's good practice to `scale` the inputs (the target is left alone). [Different techniques](https://scikit-learn.org/stable/modules/preprocessing.html) can be used, and one of the most common approaches is the z-score. Each $x_j^{(i})$ is replaced with

<div>
$$
    \frac{x_j^{(i)} - \mu_j}{\sigma_j},
$$
</div>

where $\mu_j$ is the mean of $x_j$ across all instances and $\sigma_j$ is the standard deviation. In practice, not standardizing the values can make tuning the gradient descent procedure more difficult.

## Implementation in Python

With the math written out, the implementation of the method is surprisingly simple. We just update the weights based on Equation $\eqref{eq:gradient-descent}$.

```python
import numpy as np


def scale_values(X):
    mu = np.mean(X, axis=0)
    sigma = np.std(X, axis=0)
    return (X - mu) / sigma


def linear_regression(X, y, lr, iterations, scale=True):
    # Scale the features.
    if scale:
        X = scale_values(X)

    # Prepend a column of 1s to X for the intercept term.
    n = y.shape[0]
    X_intcpt = np.hstack((np.ones((n, 1)), X))

    weights = np.zeros(X_intcpt.shape[1])
    epoch_errs = np.empty(iterations)

    # Run gradient descent for the given number of instances.
    for i in range(iterations):
        y_hat = np.matmul(X_intcpt, weights)
        epsilon = (y_hat - y).reshape(-1, 1)

        dJ = np.sum(X_intcpt * epsilon, axis=0)
        weights = weights - lr * dJ

        epoch_errs[i] = np.sum(np.square(epsilon) / n)
    return (weights, epoch_errs)
```

To test our function, we will generate a fake dataset designed for regression problems. Scikit-learn never ceases to amaze me:

```python
from sklearn.datasets import make_regression

X, y = make_regression(
    n_samples=500, n_features=2, n_informative=2, noise=2.0, random_state=42
)
```

We have a data set with 500 observations of two features that can predict $y$ pretty well. Now let's try applying our function on this data:

```python
coefs, errs = linear_regression(X, y, lr=0.0005, iterations= 500)
coefs
# array([ 0.64392818, 10.48242209, 16.19542938])
```

### Checking results

To see if these fitted coefficients are actually the ones that minimize the mean squared error, we may use the formula[^linear-regression] to compute the coefficients directly:

```python
X_scaled = scale_values(X)
X_scaled = np.hstack((np.ones((y.shape[0], 1)), X_scaled))

np.linalg.inv(X_scaled.T @ X_scaled) @ X_scaled.T @ y
# array([ 0.64392818, 10.48242209, 16.19542938])
```

Yep, exactly the same values! Actually if we plot the `errs` variable, it reaches the global minimum after less than 20 iterations.

### Using scikit-learn

Finally, we fit this multiple regression model using scikit-learn:

```python
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

scaler = StandardScaler().fit(X)
X_scaled = scaler.transform(X)
reg = LinearRegression(normalize=True)
reg.fit(X_scaled, y)
print(f"Intercept: {reg.intercept_}")
reg.coef_
# Intercept: 0.643928183919589
# array([10.48242209, 16.19542938])
```

Once again, we get the same coefficients.
