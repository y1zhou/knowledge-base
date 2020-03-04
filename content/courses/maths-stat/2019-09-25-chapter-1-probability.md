---
title: Chapter 1. Probability
author: Yi Zhou
date: '2019-09-25'
slug: chapter-1-probability
categories:
  - Statistics
tags:
  - Statistics
  - STAT 6510
linktitle: 2019 09 25 Chapter 1 Probability
summary: ~
lastmod: '2020-03-03T23:44:26-05:00'
toc: yes
type: docs
menu:
  notes:
    name: Chapter 1. Probability
    weight: 2
---

In this chapter, we introduce the concept of the probability of an event. Then we show how probabilities  can be computed in certain situations. As a preliminary, however, we need to discuss the concept of the sample space and the events of an experiment.

## Experiment, sample space and events

An experiment is the process by which an observation is made. In particular, we are interested in a `random experiment` whose outcome is not predictable with certainty. The set of *all* possible outcomes of an experiment is known as the `sample space` of the experiment, and is often denoted $S$. An `event` (denoted $E$) is a set that contains some possible outcomes of the random experiment.

By definition, any event is a subset of the sample space. For a given random experiment, its sample space is unique. Let's see some examples.

| Experiment                             | Possible outcome                 | Sample space                             | Event                                                        |
| -------------------------------------- | -------------------------------- | ---------------------------------------- | ------------------------------------------------------------ |
| Test of a certain disease on a patient | Positive of negative             | $S = \\{p, n\\}$                         | $E = \\{n\\}$                                                |
| Rolling a six-sided die                | $1, \cdots, 6$                   | $S = \\{1, 2, 3, 4, 5, 6\\}$             | Outcome $\geq 3: E = \\{4, 5, 6\\}$                          |
| Tossing two coins                      | Each coin is either head or tail | $S = \\{ (H,H), (H,T), (T,H), (T,T) \\}$ | First toss is a head: $E = \\{ (H, H), (H, T) \\}$           |
| Life time of a computer                | All non-negative real numbers    | $S = \\{X: 0 \leq X < \infty \\}$        | Survives for more than 10 hours: $E = \\{ X: 10 < X < \infty \\}$ |

For each experiment, we may define more than one event. Take the die-rolling example, we can define events like

$$
E_3 = \\{3\\} \qquad E_4 = \\{4\\} \qquad E_5 = \\{5\\} \qquad E_6 = \\{6\\}
$$

If we observed the event $E = \\{4, 5, 6\\}$, it means we observed one of the three events $E_4$, $E_5$ or $E_6$. We say $E$ can be decomposed into $E_4$, $E_5$, and $E_6$. If an event can be further decomposed, it is called a `compound event`. Otherwise, it's called a `simple event`. Each simple event contains one and only one outcome.

Finally, events with no outcome is called the `null event`, and is denoted $\emptyset$. For example, an event of the outcome is greater than 7 in the die-rolling example.

