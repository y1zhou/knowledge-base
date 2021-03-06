---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Basic Concepts"
slug: mathematical-statistics-basic-concepts
summary: "Introducing the concept of the probability of an event. Also covers set operations and the sample-point method."
categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
date: 2019-09-25T11:05:06-05:00
toc: true # Show table of contents? true/false
featured: true
weight: 11

header_image:
  unsplash_id: "mbXEkW5ZyBQ"
  caption: "Macro of a Intel motherboard"
  
---

In this chapter, we introduce the concept of the probability of an event. Then we show how probabilities can be computed in certain situations. As a preliminary, however, we need to discuss the concept of the sample space and the events of an experiment.

## Experiment, sample space and events

An experiment is the process by which an observation is made. In particular, we are interested in a `random experiment` whose outcome is not predictable with certainty. The set of _all_ possible outcomes of an experiment is known as the `sample space` of the experiment, and is often denoted $S$. An `event` (denoted $E$) is a set that contains some possible outcomes of the random experiment.

By definition, any event is a subset of the sample space. For a given random experiment, its sample space is unique. Let's see some examples.

| Experiment                             | Possible outcome                 | Sample space                             | Event                                                             |
| -------------------------------------- | -------------------------------- | ---------------------------------------- | ----------------------------------------------------------------- |
| Test of a certain disease on a patient | Positive of negative             | $S = \\{p, n\\}$                         | $E = \\{n\\}$                                                     |
| Rolling a six-sided die                | $1, \cdots, 6$                   | $S = \\{1, 2, 3, 4, 5, 6\\}$             | Outcome $\geq 3: E = \\{4, 5, 6\\}$                               |
| Tossing two coins                      | Each coin is either head or tail | $S = \\{ (H,H), (H,T), (T,H), (T,T) \\}$ | First toss is a head: $E = \\{ (H, H), (H, T) \\}$                |
| Life time of a computer                | All non-negative real numbers    | $S = \\{X: 0 \leq X < \infty \\}$        | Survives for more than 10 hours: $E = \\{ X: 10 < X < \infty \\}$ |

For each experiment, we may define more than one event. Take the die-rolling example, we can define events like

$$
\begin{aligned}
	&E_3 = \\{3\\} &E_4 = \\{4\\} \\\\
  &E_5 = \\{5\\} &E_6 = \\{6\\}
\end{aligned}
$$

If we observed the event $E = \\{4, 5, 6\\}$, it means we observed one of the three events $E_4$, $E_5$ or $E_6$. We say $E$ can be **decomposed** into $E_4$, $E_5$, and $E_6$. If an event can be further decomposed, it is called a `compound event`. Otherwise, it's called a `simple event`. Each simple event contains one and only one outcome.

Finally, events with no outcome is called the `null event`, and is denoted $\emptyset$. For example, an event of the outcome is greater than 7 in the die-rolling example.

## Set operations

Suppose we have a sample space $S$ and two events $E$ and $F$.

### Subset

If all of the outcomes in $E$ are also in $F$, then we say that $E$ is contained in $F$, or $E$ is a `subset` of $F$. We write it as $E \subset F$. Subsets have several properties:

1. Any event is a subset of the sample space: $E \subset S$.
2. Any event is a subset of itself: $E \subset E$.
3. If $E \subset F$ and $F \subset E$, then $E = F$.
4. $\emptyset \subset E$, $\emptyset \subset S$.
5. $E \subset F, F \subset G \Rightarrow E \subset G$.

### Union

We denote $E \cup F$ as the `union` of the two events. It is a new event which consists of all outcomes in $E$, and all the outcomes in $F$. In other words, $E \cup F = $ {either in $E$ or in $F$}. The union operation also has a few properties:

1. $E \cup F \subset S$.
2. $E \cup E = E$.
3. $E \cup S = S, E \cup \emptyset = E$.

### Intersect

We denote $E \cap F$, or $EF$ for short, the `intersection` between $E$ and $F$. $E \cap F$ consists of all outcomes that are both in $E$ and in $F$. Its properties are straightforward:

1. $E \cap S = E$.
2. $E \cap E = E$.
3. $E \cap \emptyset = \emptyset$.

We say $E$ and $F$ are `disjoint` or `mutually exclusive` if $E \cap F = \emptyset$. Any event is disjoint with the null event.

### Complement

Finally, for any event $E$, we define a new event $E^C$, referred to as the `complement` of $E$. $E^C$ consists of all outcomes in the sample space $S$ that are not in $E$. Its properties include

1. $S^C = \emptyset$ and $\emptyset^C = S$.
2. $E \cup E^C = S$.
3. The event $E$ is always disjoint with its complement: $E \cap E^C = \emptyset$.
4. $\left( E^C \right)^C = E$.

### Example of set operations

Consider rolling two six-sided dice. Let

$$
\begin{aligned}
	E_1 &= \\{\text{first roll is } 3 \\} \\\\
	E_2 &= \\{\text{sum of two rolls is } 7 \\} \\\\
	E_3 &= \\{\text{second roll} - \text{first roll} < 4 \\}
\end{aligned}
$$

and we want to find (1) $E_1 \cup E_2$, (2) $E_1 \cap E_2$, and (3) $E_3^C$.

The sample space is

$$
S = \\{\underbrace{ (1, 1), \cdots, (6, 6)}_\text{36 outcomes} \\}
$$

and the events are

$$
\begin{aligned}
	E_1 &= \\{ (3, 1), (3, 2), \cdots, (3, 6) \\} \\\\
	E_2 &= \\{ (1, 6), (2, 5), (3, 4), (4, 3), (5, 2), (6, 1) \\}
\end{aligned}
$$

We skip $E_3$ for now as it's more complicated. We can find that $E_1$ and $E_2$ only have one element in common, $(3, 4)$, so $E_1 \cap E_2 = \\{ (3, 4) \\}$, and $E_1 \cup E_2 = $ {$(3, 1), (3, 2), (3, 3), (3, 5), (3, 6), E_2$}.

Finally, we can express the complement of $E_3$ as $E_3^C = \\{ \text{second} - \text{first} \geq 4 \\}$, which is $\\{ S - F = 4 \\} \cup \\{ S - F = 5 \\}$, and we can write down all the possibilities as $\\{ (1, 5), (1, 6), (2, 6) \\}$.

A graphical representation that is useful for illustrating logical relations among events is the [Venn diagram](https://en.wikipedia.org/wiki/Venn_diagram).

### Laws of set operations

The operations of sets can be applied to more than two events, and they follow certain rules similar to the rules of algebra. All the following rules can be verified by Venn diagrams.

- Commutative laws: $E \cup F = F \cup E$, $E \cap F = F \cap E$.
- Associative laws: $(E \cup F) \cup G = E \cup (F \cup G)$, $(E \cap F) \cap G = E \cap (F \cap G)$.
- Distributive laws: $(E \cup F) \cap G = ( E \cap G ) \cup ( F \cap G )$, $(E \cap F) \cup G = (E \cup G) \cap (F \cup G)$.

In addition, there is a law that connects all three operations (union, intersection and complement) together, `DeMorgan's law`:

$$
\begin{aligned}
	(E \cup F)^C &= E^C \cap F^C \\\\
	(E \cap F)^C &= E^C \cup F^C
\end{aligned}
$$

DeMorgan's law can be extended to more than two events. Let $\bigcup_{i=1}^n E_i$ denote the union of events $E_1$ to $E_n$, and $\bigcap_{i=1}^n E_i$ their intersection,

$$
\begin{aligned}
	\left( \bigcup_{i=1}^n{E_i} \right)^C &= \bigcap_{i=1}^n{E_i^C} \\\\
	\left( \bigcap_{i=1}^n{E_i} \right)^C &= \bigcup_{i=1}^n{E_i^C}
\end{aligned}
$$

## Probability of events

One way of defining the probability of an event is in terms of its relative frequency. Suppose we have a random experiment with sample space $S$, and we want to assign some number $P(E)$ to represent the probability of event $E$. We may repeat this random experiment many times. Let $n(E)$ be the number of times in the first $n$ repetitions of the experiment that the event $E$ occurs. The probability of the event is defined as

$$
P(E) = \lim_{n \rightarrow \infty}\frac{n(E)}{n}
$$

There are a few drawbacks to this method:

1. It requires $S$ to be countable.
2. We need to assume that the limit exists, and is a positive number.
3. Sometimes our random experiments are limited and we can't repeat it many times, or the experiment may not even be observable.

To overcome such drawbacks, modern mathematics used an axiom system to define the probability of an event.

### Axioms of probability

For sample space $S$ and event $E$, we define three axioms.

**Axiom 1.** $0 \leq P(E) \leq 1$.

**Axiom 2.** $P(S) = 1$.

**Axiom 3.** For any sequence of mutually exclusive events $E_1, E_2, \cdots$,

$$
P\left(\bigcup_{i=1}^\infty E_i \right) = \sum_{i=1}^\infty{P(E_i)}
$$

where $E_i \cap E_j = \emptyset$ for any $i$ and $j$ where $i \neq j$ (mutually exclusive). More formally, we can say $P$ to be $\sigma$-additive.

The definition through axioms is mathematically rigorous, flexible, and can be developed into an axiomatic system. We'll show the flexibility through an example. Suppose our experiment is tossing a coin. If we believe it is a fair coin, we have

$$
S = \\{ H, T \\}, \quad P(\\{H\\}) = P(\\{T\\})
$$

then using Axioms $2$ and $3$ above, we can derive

$$
P\left(\bigcup_{i=1}^2 E_i \right) = \sum_{i=1}^2{P(E_i)} = P(S) = 2P(\\{H\\}) = 2P(\\{T\\})
$$

so $P(\\{H\\}) = P(\\{T\\}) = 0.5$. If we believe the coin is biased and $P(\\{H\\}) = 2P(\\{T\\})$, then

$$
P(\\{H\\}) + P(\\{T\\}) = P(\\{H\\} \cup \\{T\\}) = P(S) = 1
$$

By combining the two equations, $3P(\\{T\\}) = 1 \Rightarrow P(\\{T\\}) = \frac{1}{3}$, $P(\\{H\\}) = \frac{2}{3}$.

In this example, we didn't use any information of the observations or frequencies. We are assigning probabilities according to our belief so long as this assignment satisfies the three axioms. Based on the axioms, we can prove some simple propositions of probability.

### Propositions

**Proposition 1:** $P(E^C) = 1 - P(E)$. The proof is given as follows.

$$
\begin{gather*}
	E^C \cap E = \emptyset,\\, E^C \cup E = S \\\\
	P(E^C \cup E) = P(E^C) + P(E) = 1 \\\\
  \Rightarrow P(E^C) = 1 - P(E)
\end{gather*}
$$

**Proposition 2:** If $E \subset F$, then $P(E) \leq P(F)$. To prove this, note that $E$ and $E^C \cap F$ are mutually exclusive.

$$
\begin{aligned}
	E \cap (F \cap E^C) &= E \cap (E^C \cap F) = (E \cap E^C) \cap F = \emptyset \cap F = \emptyset \\\\
  F &= E \cup (F \cap E^C) \\\\
  P(F) &= P(E) + \underbrace{P(F \cap E^C)}_{\geq 0} \geq P(E)
\end{aligned}
$$

![Venn diagram of Proposition 2.](1-proposition-2.svg)

**Proposition 3:** $P(E \cup F) = P(E) + P(F) - P(EF)$. This proposition can be easily proved using a Venn diagram. Let $I$, $II$ and $III$ denote $E \cap (F^C)$, $E \cap F$ and $F \cap E^C$, respectively.

$$
\begin{aligned}
	P(E) &= P(I) + P(II) \\\\
  P(F) &= P(II) + P(III) \\\\
  P(E \cup F) &= P(I) + P(II) + P(III) \\\\
  &= P(E) + P(F) - P(II) \\\\
  &= P(E) + P(F) - P(EF)
\end{aligned}
$$

![Venn diagram of Proposition 3.](1-proposition-3.svg)

Now let's apply these propositions to the example below. A student is applying for two jobs. Suppose she'll get an offer from company A with probability $0.3$, and an offer from company B with probability $0.4$, and with probability $0.3$ she gets both offers. What is the probability that she gets neither offer?

$$
\begin{aligned}
	S &= \\{(S, S), (S, F), (F, S), (F, F)\\} \\\\
  E &= \\{\text{get offer from A}\\} = \\{(S, S), (S, F)\\} \\\\
  F &= \\{\text{get offer from B}\\} = \\{(S, S), (F, S)\\} \\\\
  G &= \\{\text{get two offers}\\} = \\{(S, S)\\} \\\\
  K &= \\{\text{get no offers}\\} = \\{(F, F)\\}
\end{aligned}
$$

We have $G = E \cap F$ and $K = E^C \cap F^C = (E \cup F)^C$. Knowing that $P(E) = 0.5, P(F) = 0.4$ and $P(G) = P(E \cap F) = 0.3$,

$$
\begin{aligned}
	P(K) &= P\left( (E \cup F)^C \right) \\\\
	&=  1 - P(E \cup F) = 1 - (P(E) + P(F) - P(EF)) = 0.4
\end{aligned}
$$

## The sample-point method

As somewhat shown in [the example above](#propositions "propositions example"), for an experiment with finite or countable number of outcomes, we can calculate the probability of an event through the so called `sample-point method`. The procedure is

1. Define the experiment, sample space and simple events (outcomes).
2. Assign reasonable probabilities to each simple event.
3. Define the event of interest as a collection of simple events.
4. Calculate the probability of the event by summing the probabilities of the simple events in the event.

The main idea of the sample-point method is based on Axiom $3$.

### Coin flip example

A fair coin is tossed three times. Find the probability that _exactly_ two of the three tosses are heads.

We'll follow the procedure in the sample-point method. The experiment is "tossing the coin three times". The sample space is

$$
S = \\{\underbrace{(H, H, H), (T, H, H), \cdots, (T, T, T)}_8\\}
$$

where each of the $8$ outcomes can be considered as a simple event $E_1, \cdots, E_8$. Since we consider it as a fair coin,

$$
P(E_1) = P(E_2) = \cdots = P(E_8) = \frac{1}{8}
$$

Our event of interest, $E$, is defined as {2 heads and 1 tail}, so

$$
\begin{aligned}
	E &= \\{(T, H, H), (H, T, H), (H, H, T)\\} \\\\
  &= \\{(T, H, H)\\} \cup \\{(H, T, H)\\} \cup \\{(H, H, T)\\}
\end{aligned}
$$

The final step is calculating the probability of $E$

$$
P(E) = P(F_1 \cup F_2 \cup F_3) = \sum_{i=1}^3P(F_i) = \frac{3}{8}
$$

Note that in this case (and in many other experiments), all the outcomes in the sample space are equally likely to occur. For such experiments, we can simplify the sample-point method as

$$
P(E) = \frac{\text{\# of outcomes in }E}{\text{\# of outcomes in }S}
$$

### Powerball example

The Powerball is one of the largest lottery games in the US. The system works like this:

1. $5$ numbered white balls are drawn out of $69$ balls without replacement.
2. $1$ numbered red ball is drawn out of $26$ balls.

You win the Powerball if you chose exactly those $5+1$ balls, and the order of the white balls doesn't matter. What is the probability to win a Powerball?

It's reasonable to assume each outcome will be equally likely to occur. Each outcome is a set of $6$ numbers satisfying the above rules.

$$
\begin{aligned}
    E &= \{ \text{You win the PB} \} \\\\
    &= \{\text{You choose exactly the lucky number}\} \\\\
    |E| &= 1 \\\\
    P(E) &= \frac{|E|}{|S|} = \frac{1}{|S|} \\\\
    |S| &= \text{\# ways to draw 1 red ball out of 26 } (26) \\\\
    &\quad \times \text{\# ways to draw 5 white balls out of 69}
\end{aligned}
$$

If we draw the white balls one by one, we have $69 \times 68 \times 67 \times 66 \times 65$ ways (ordered outcomes). For each set of $5$ numbers, we have $5 \times 4 \times 3 \times 2$ ways of arranging them. So the number of ways to draw $5$ white balls of $69$ is

$$
\frac{69 \times 68 \times 67 \times 66 \times 65}{5 \times 4 \times 3 \times 2}
$$

Formally, this is "choose $k$ from $n$", which can be written as $\binom{n}{k}$ and

$$
\binom{n}{k} = \frac{n(n-1)\cdots(n-k+1)}{k(k-1)(k-2)\cdots} = \frac{n!}{(n-k)!k!}
$$

Now we have

$$
\begin{gather*}
    |S| = 26 \times \binom{69}{5} = 26 \times 11,238,513 \approx 292M \\\\
    P(E) = \frac{1}{|S|} \approx 3.42 \times 10^{-9}
\end{gather*}
$$
