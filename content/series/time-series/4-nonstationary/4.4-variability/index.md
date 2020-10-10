---
title: "Variability of Nonstationary Time Series"
date: 2020-10-09T11:30:59-04:00
summary: "" # appears in list of posts
categories: ["Time Series"] # main category; shown in post metadata
tags: ["Time Series", "Visualization", "R"] # list of related tags

slug: "time-series-stationarity-variability"
toc: true # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 44 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "Beer collection." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "nhX8QhXMBkM" # Unsplash ID of the picture
---

In most of the chapter we've been talking about the mean trend of a time series and various methods of detecting and removing it. However, there's another type of nonstationarity where the variability of a time series change with the mean level over time.

Recall that in linear regression we apply a log-transformation on $y$ when the variance of the residuals increases as $\hat{y}$ increases. Similar techniques can be applied for time series data. Here we're talking about global changes -- if the variance changes locally, then a transformation wouldn't help and different methods would be needed to stabilize the variance.

## Log transformation

We'll use the `AirPassengers` dataset in R to show how the variability of a time series could be stabilized.

{{< figure src="air_passengers.png" caption="Monthly totals of international airline passenger numbers from 1949 to 1960." numbered="true" >}}

We can see clearly from the time series plot that the mean and variance both increase over time, and there's a seasonal pattern involved. The increasing mean trend can be handled by differencing, and the seasonal pattern will be discussed in the next chapter.

For stabilizing the variance, we can log-transform the dataset. As we can see below, the mean still increases, but the variance becomes stable across the years.

{{< figure src="log_air_passengers.png" caption="Log-transformed Monthly totals of international airline passenger numbers from 1949 to 1960." numbered="true" >}}
