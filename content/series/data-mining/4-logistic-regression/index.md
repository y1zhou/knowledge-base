---
title: "Logistic Regression"
date: 2021-02-12T15:55:00-04:00
summary: "" # appears in list of posts
categories: ["Data Mining"] # main category; shown in post metadata
tags: ["Data Mining", "Regression"] # list of related tags

slug: "data-mining-logistic-regression"
toc: false # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 40 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Superman logo." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "6RcT6zEmm9A" # Unsplash ID of the picture
---

In linear regression, the function learned is used to estimate the value of the target $y$ using values of input $x$. While it could be used for classification purposes by setting the target value to a distinct constant for each class, it's a poor choice for this task. The target attribute takes on a finite number of values, yet the linear model produces a continuous range.

For classification tasks, `logistic regression` is a better choice. As with linear regression, logistic regression learns weights for a linear equation. However, instead of using the equation to predict a target attribute's value, it separates instances into classes.

## Maths of logistic regression

### Hypothesis function

Simple logistic regression works on two classes, i.e. $y=0$ or $y=1$. Rather than using `$h_\theta (x) = \sum_{j=0}^k \theta_j x_j$` as the hypothesis, we use a hypothesis based on the `sigmoid` or `logistic function`, which is continuous and bounded by 0 and 1:

<div>
$$
    f(z) = \frac{1}{1 + e^{-z}}
$$
</div>

Note that $f(z)$ increases as $z$ increases, and $f(0) = 0.5$. Multiplying $z$ by a constant affects how steep the curve is, while adding a constant to $z$ shifts the curve to the left or right. The hypothesis function used in logistic regression is:

<div>
$$
    h_\theta (x) = \frac{1}{1 + e^{-\sum_{j=0}^k \theta_j x_j}}
$$
</div>

The values of $h_\theta (x)$ are probabilities. Once the proper weights $\theta$ are learned, `$\sum_{j=0}^k \theta_j x_j = 0$` defines the `decision boundary` for the classifier. Any instance $x$ where `$\sum_{j=0}^k \theta_j x_j \geq 0$` falls on one side of the boundary, and `$\sum_{j=0}^k \theta_j x_j < 0$` falls on the other.

### Cost function

In this context, using $h_\theta$ in the cost function used for linear regression results in a function that's not convex. To ensure a convex function, the following function is used:

<div>
$$
    J(\theta) = \sum_{i=1}^m -y^{(i)} \log \left(h_\theta (x^{(i)}) \right) - (1 - y^{(i)}) \log \left(1 - h_\theta (x^{(i)}) \right)
$$
</div>

This is actually the combination of two functions. Given one instance $(x, y)$,

<div>
$$
    J(\theta) = \begin{cases}
        -\log(h_\theta(x)), & y=1 \\
        -\log(1 - h_\theta (x)), & y=0 
    \end{cases}
$$
</div>

The idea is that if $y=1$, then the cost function approaches 0 as $h_\theta(x)$ approaches 1, and as $h_\theta (x)$ approaches 0, the cost function increases. The case is similar when $y=0$. In other words, if the actual and predicted values are the same, then the error is 0. However, as the predicted value moves away from the actual one, the error grows. Note how the cost is always non-negative.

### Gradient descent

Given the hypothesis and the cost function, we can work out that the update for each weight $\theta_j$ in gradient descent is the same as it is in linear regression:

<div>
$$
    \theta_j' = \theta_j - \alpha \left( \sum_{i=1}^n \left[ (h_\theta (x^{(i)}) - y^{(i)}) x_j^{(i)} \right] \right)
$$
</div>

## Implementation in Python

The code for logistic regression is really similar to [that of linear regression's]({{< ref "/series/data-mining/3-gradient-descent-linear-regression/index.md#implementation-in-python" >}}).

```python
import numpy as np

def logistic_regression(X, y, lr, iterations):
    # Prepend a column of 1s to X for the intercept term.
    n = y.shape[0]
    X_intcpt = np.hstack((np.ones((n, 1)), X))

    weights = np.zeros(X_intcpt.shape[1])
    epoch_accuracy = np.empty(iterations)

    # Run gradient descent for the given number of instances.
    for i in range(iterations):
        y_hat = np.reciprocal(1 + np.exp(-1 * np.matmul(X_intcpt, weights)))
        epsilon = (y_hat - y).reshape(-1, 1)

        dJ = np.sum(X_intcpt * epsilon, axis=0)
        weights = weights - lr * dJ

        predicted_class = np.where(y_hat >= 0.5, 1, 0)
        epoch_accuracy[i] = np.sum(predicted_class == y) / n

    return (weights, epoch_accuracy)
```

To test this function, we generate a classification dataset using scikit-learn, and compare the prediction results:

```python
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

X, y = make_classification(
    n_samples=500, n_features=2, n_informative=2, n_redundant=0, random_state=42
)

lr = 2e-5
iterations = 30000

coefs, accuracies = logistic_regression(X, y, lr, iterations)
print(f"Accuracy: {accuracies[-1]}")
coefs
# Accuracy: 0.892
# array([-0.17091739,  2.58437914,  0.25970821])

# scikit-learn results
reg = LogisticRegression()
reg.fit(X, y)
print(f"Accuracy: {reg.score(X, y)}")
print(f"Intercept: {reg.intercept_}")
# reg.coef_
# Accuracy: 0.89
# Intercept: [-0.16069779]
# array([[2.46293031, 0.24200195]])
```

The values are pretty close, although our method takes much longer to run.
