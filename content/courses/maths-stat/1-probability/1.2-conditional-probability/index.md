---
# Documentation: https://sourcethemes.com/academic/docs/managing-content/

title: "Conditional Probability"
slug: mathematical-statistics-conditional-probability
summary: "Introducing conditional probability and independence of events. Bayes' rule comes in as well."
date: 2019-09-26T11:51:56-05:00

categories:
  - Mathematical Statistics
tags:
  - Mathematical Statistics
  - Statistics
  - Bayes

toc: true  # Show table of contents? true/false
type: docs  # Do not modify.
weight: 20

menu:
  maths-stat:
    name: Conditional Probability
    parent: Probability
    weight: 20
---

## Conditional probability and independence of events

The probability of an event will sometimes depend upon whether we know that other events have occurred. This is easier to explain with an example.

### Conditional probability

Suppose we roll two fair six-sided dice. What is the probability that the sum of the two dice is 8? Using the procedure above, we can easily get

$$
P(E) = \frac{|E|}{|S|} = \frac{|\{(2, 6), (3, 5), (4, 4), (5, 3), (6, 2)\}|}{6 \times 6} = \frac{5}{36}
$$

What if we know the first roll is 3? That would be another event:

$$
\begin{aligned}
    F &= \\{\text{first roll is } 3\\} \\\\
    E' &= \\{ \text{sum is } 8 \text{ given } F \\}
\end{aligned}
$$

Given the first die is 3, requiring the sum to be 8 is equivalent to requiring the second roll to be 5. So the probability is

$$
P(E') = \frac{|E'|}{|S'|} = \frac{|\{(3, 5)\}|}{|\{(3, 1), \cdots, (3, 6)\}|} = \frac{1}{6} \neq \frac{5}{36}
$$

Formally speaking, let $E$ be the event that sum is 8, and $F$ be the event that the first roll is 3. The `conditional probability` of $E$ given $F$ is denoted $P(E \mid F)$, and $P(E)$ is the unconditional probability of $E$. If $P(F) > 0$, then

$$
\begin{equation}P(E \mid F) = \frac{P(EF)}{P(F)} \label{eq:conditional-prob} \end{equation}
$$

To understand this, keep in mind that any event $E$ can be decomposed into $(EF) \cup (EF^C)$. From Eq. $\eqref{eq:conditional-prob}$, we can derive

$$
\begin{equation} P(EF) = P(E \mid F)P(F) \end{equation}
$$

Now we can revisit the example above:

$$
\begin{gather*}
    \text{Let } E = \\{\text{sum} = 8\\}, \quad F = \\{ \text{first}\\} = 3. \\\\
    P(E \mid F) = \frac{P(EF)}{P(F)} \\\\
    P(EF) = \frac{|\\{(3, 5)|\\}}{|S|} = \frac{1}{36} \\\\
    P(F) = \frac{|\\{ (3, 1), \cdots, (3, 6) \\}|}{|S|} = \frac{6}{36}
\end{gather*}
$$

The conditional probability can also be generalized to more than two events using the `multiplication rule`:

$$
\begin{equation}\begin{split}
    &\quad P(E_1)P(E_2 \mid E_1)P(E_3 \mid E_1 E_2) \cdots P(E_n \mid E_1 E_2 \cdots E_{n-1}) \\\\
    &= P(E_1) \frac{P(E_1E_2)}{P(E_1)} \frac{P(E_3E_1E_2)}{P(E_1E_2)} \cdots \frac{P(E_nE_1E_2 \cdots E_{n-1})}{P(E_1E_2 \cdots E_{n-1})} \\\\
    &= P(E_1 E_2 \cdots E_n)
\end{split} \end{equation}
$$

### Cards example

Suppose we have a deck of 52 cards, and we randomly divided them into 4 piles of 13 cards. Compute the probability that each pile has exactly 1 ace.

There are four houses in a deck of cards: Hearts, Diamonds, Clubs and Spades. We can define events $E_i, i = 1, 2, 3, 4$ as follows

$$
\begin{aligned}
    E_1 &= \\{ \text{Ace of Hearts is in any pile} \\} \\\\
    E_2 &= \\{ \text{Ace of Hearts and Diamonds are in different piles} \\} \\\\
    E_3 &= \\{ \text{Ace of Hearts, Diamonds and Clubs are in different piles} \\} \\\\
    E_4 &= E = \\{ \text{All four aces are in different piles} \\}
\end{aligned}
$$

The desired probability is $P(E_4)$.

$$
\begin{aligned}
    P(E_4) &= P(E_1 E_2 E_3 E_4) \\\\
    &= P(E_1) P(E_2 \mid E_1) P(E_3 \mid E_1E_2) P(E_4 \mid E_1E_2E_3) \\\\
    &= P(E_1) P(E_2 \mid E_1) P(E_3 \mid E_2) P(E_4 \mid E_3)
\end{aligned}
$$

$P(E_1) = 1$ because the ace of Hearts is always going to be in a pile. For $P(E_2 \mid E_1)$, we can calculate the probability of the two cards being in the same pile. As the remaining 12 cards are equally likely drawn from the deck of 51 cards,

$$
P(E_2 \mid E_1) = 1 - \frac{12}{51} = \frac{39}{51}
$$

Given $E_1E_2$, the ace of clubs can't be any of the 24 cards in the two piles with the two aces.

$$
P(E_3 \mid E_2) = 1 - \frac{12 + 12}{50} = \frac{26}{50}
$$

and finally we have

$$
P(E_4 \mid E_3) = 1 - \frac{12 \times 3}{49} = \frac{13}{49}
$$

So $P(E_4) = 1 \cdot \frac{39}{51} \cdot \frac{26}{50} \cdot \frac{13}{49} \approx 0.105$.

### Independence of events

In general, the conditional probability $P(E \mid F) \neq P(E)$. We say $E$ and $F$ are independent when $P(E \mid F) = P(E)$. When $E$ is independent of $F$, we also have

$$
P(EF) = P(E)P(F)
$$

Two events $E$ and $F$ are `independent` if any of the following holds:

$$
\begin{cases}
    P(E \mid F) = P(E) \\\\
    P(F \mid E) = P(F) \\\\
    P(EF) = P(E)P(F)
\end{cases}
$$

Otherwise we say $E$ and $F$ are dependent. Independence is denoted by $E \perp F$.

**Proposition:** if $E$ and $F$ are independent, then so are $E$ and $F^C$.

**Proof:** we need to find $P(EF^C) = P(E)P(F^C)$.

$$
\begin{equation} \label{eq:event-independence}
\begin{split}
    E = E \cap S &= E \cap (F \cup F^C) \\\\
    &= (E \cap F) \cup (E \cap F^C) \\\\
    P(E) &= P(EF) + P(EF^C) \\\\
    &= P(E)P(F) + P(EF^C) \\\\
    P(E) - P(E)P(F) &= P(EF^C) \\\\
    P(E)P(F^C) &= P(EF^C)
\end{split}
\end{equation}
$$

From this we have

$$
E \perp F \Rightarrow E \perp F^C \Rightarrow E^C \perp F \Rightarrow E^C \perp F^C
$$

We'll finish up this part with another example. Suppose we have a circuit with $n$ switches in parallel. The probability that the $i^{th}$ component can work is $P_i$, $i = 1, \cdots, n$. What is the probability that the system functions?

Denote $E_i =$ {the $i^{th}$ component functions}, $P(E_i) = P_i$, and $E_i \perp E_j \, \forall i \neq j$.

$$
\begin{gather*}
    E = \\{\text{system functions}\\} = \\{\text{observe at least one } E_i\\} \\\\
    E^C = \\{\text{none of the components can work}\\} = \\{\text{didn't observe any } E_i\\} = \bigcap_{i=1}^n{E_i^C} \\\\
    P\left(E^C\right) = P\left(\bigcap_{i=1}^n{E_i^C}\right) = \prod_{i=1}^n{P\left(E_i^C\right)} \\\\
    P(E) = 1 - P(E^C) = 1 - \prod_{i=1}^n{P\left(E_i^C\right)}
\end{gather*}
$$

## Law of total probability and Bayes' rule

In Proof $\eqref{eq:event-independence}$, we showed a trick to represent the probability of an event

$$
P(E) = P(EF) + P(EF^C)
$$

because $F \cup F^C = S$ and $FF^C = \emptyset$. Now let's consider a generalization of this. For some positive integer $k$, let the sets $E_1, \cdots, E_k$ be such that $\bigcup_{i=1}^k {E_i} = S$ and $E_i \cap E_j = \emptyset \quad \forall i \neq j$. Then the collection of sets $E_1, \cdots, E_k$ is called a `partition` of $S$. $F$ and $F^C$ is a partition of $S$ with $k=2$.

### Law of total probability

The theorem states that given $F_1, \cdots, F_k$ as a partition of $S$, such that $P(F_i) > 0$ for $i = 1, \cdots, k$, for any event $E$, we have

$$
P(E) = \sum_{i=1}^k{P(E \mid F_i)P(F_i)}
$$

and the proof is given as follows.

$$
\begin{aligned}
    E &= E \cap S = E \cap \left( \bigcup_{i=1}^k{F_i} \right) \\\\
    &= (E \cap F_1) \cup (E \cap F_2) \cup\cdots\cup (E \cap F_k) \\\\
    P(E) &= P\left( \bigcup_{i=1}^k{EF_i} \right) \\\\
    &= \sum_{i=1}^k{P(EF_i)} \quad \text{because } EF_i \text{ are pairwise disjoint} \\\\
    &= \sum_{i=1}^k{P(E \mid F_i)P(F_i)}
\end{aligned}
$$

**Example:** In a driving behavior survey, $60\%$ are sedan drivers, $30\%$ are SUV drivers, and $10\%$ are other drivers. $40\%, 65\%$ and $55\%$ of sedan, SUV and other drivers have received a citation within the past 3 years. Suppose each driver can own one type of car, what is the probability that a randomly selected driver received a citation within 3 years?

$D_1$, $D_2$ and $D_3$ are the events that a random driver drives a sedan, SUV or other car, respectively. Let $Y$ be the event that the driver received a citation within 3 years. We have

$$
P(Y) = \sum_{i=1}^3{P(Y \mid D_i)P(D_i)}
$$

because $D_1, D_2$ and $D_3$ form a partition of $S$. 

$$
\begin{gather*}    P(D_1) = 0.6 && P(D_2) = 0.3 && P(D_3) = 0.1 \\\\    P(Y \mid D_1) = 0.4 && P(Y \mid D_2) = 0.65 && P(Y \mid D_3) = 0.55\end{gather*}
$$

Therefore $P(Y) = 0.49$.

### Bayes' rule

Using the law of total probability, we can derive a simple but very useful result known as the Bayes' rule. Assume that $F_1, \cdots, F_k$ is a partition of $S$ such that $P(F_i) > 0$ for $i = 1, \cdots, k$, then

$$
P(F_i \mid E) = \frac{P(E \mid F_i)P(F_i)}{\sum_{i=1}^k{P(E \mid F_i)P(F_i)}}
$$

By definition of conditional probability, $P(E \mid F_i)P(F_i) = P(EF_i)$. By law of total probability, we have $\sum_{i=1}^k{P(E \mid F_i)P(F_i)} = P(E)$. So

$$
\frac{P(E \mid F_i) P(F_i)}{\sum_{i=1}^k{P(E \mid F_i) P(F_i)}} = \frac{P(EF_i)}{P(E)} = P(F_i \mid E)
$$

If we only have two events, another form of Bayes' rule is

$$
P(E \mid F) = \frac{P(F \mid E)P(E)}{P(F)}
$$

if $P(F) > 0$ and $P(E) > 0$. 

With Bayes' rule, we can solve some very unintuitive conditional probability problems such as the [Monty Hall problem](https://brilliant.org/wiki/monty-hall-problem/). Here we illustrate its power with another example. A biomarker was developed to detect a certain kind of gene defect. When this test is applied to a person with this gene defect, it has a probability of $0.9$ to give a positive result. If this test is applied to a person without the defect, there's a probability of $0.05$ for the biomarker to give a false positive result. We know $1\%$ of the total population have this defect. When we apply this to a random person, what are the probabilities of

1. the test result is negative,
2. the person has the defect given the test result is positive, and
3. the person doesn't have this defect given the test result is negative.

$P$ and $N$ are the events of positive and negative results. $W$ and $O$ are the events of with and without this gene defect. We want to find $P(N)$, $P(W \mid P)$ and $P(O \mid N)$ knowing that

$$
\begin{cases}
    P(P \mid W) = 0.9 \\\\
    P(P \mid O) = 0.05 \\\\
    P(W) = 0.01
\end{cases}
$$

$$
\begin{aligned}
    P(N) &= 1 - P(P) \\\\
    &= 1 - (P(P \mid W)P(W) + P(P \mid O)P(O))) \\\\
    &= 1 - (0.9 \times 0.01 + 0.05 \times (1 - 0.01)) \\\\
    &= 1 - 0.0585 = 0.9415 \\\\
    P(W \mid P) &= \frac{P(P \mid W)P(W)}{P(P)} \\\\
    &= \frac{0.9 \times 0.01}{0.0585} = 0.1538 \\\\
    P(O \mid N) &= \frac{P(N \mid O)P(O)}{P(N)} \\\\
    &= \frac{(1 - 0.05) \times (1 - 0.01)}{0.9415} = 0.9989
\end{aligned}
$$

Next, we move on to discuss [discrete random variables]({{< ref "/courses/maths-stat/2-discrete-random-variables/2.1-definition/index.md" >}}).