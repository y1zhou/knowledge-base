---
title: "Text Classification with Naive Bayes and NLTK"
date: 2021-02-09T15:55:00-04:00
summary: "" # appears in list of posts
categories: ["Data Mining"] # main category; shown in post metadata
tags: ["Data Mining", "Classification"] # list of related tags

slug: "data-mining-text-classification-naive-bayes-nltk"
toc: false # table of contents button in post

# featured posts are shown on the homepage
featured: false
draft: false

weight: 20 # smaller values are listed first

# full-width featured image
# To use, add an image named `featured.jpg/png` to your page's folder, or
# fill the unsplash_id and the photo will be automatically retrieved.
header_image:
    caption: "A magazine stand in Kiev, Ukraine." # Give credits here, or whatever captions you want to add (support markdown)
    unsplash_id: "AM-Tkk_dkNU" # Unsplash ID of the picture
---

In the last post we talked about the theoretical side of naive Bayes in text classification. Here we will implement the model in Python, both from scratch and utilizing existing packages.

The corpus we use is a 26-line poem by T.S. Eliot. In each line a dummy string "ZZZ" or "XXX' has been inserted, representing the class of the line ("ZZZ" for class 0 and XXX for class 1).

```python
corpus = [
    "And indeed there will be time ZZZ",
    "For the yellow smoke that slides along the street XXX",
    "Rubbing its back upon the window-panes ZZZ",
    "There will be time, there will be time ZZZ",
    "To prepare a face to meet the faces that you meet XXX",
    "There will be time to murder and create ZZZ",
    "And time for all the works and days of hands ZZZ",
    "That lift and drop a question on your plate ZZZ",
    "Time for you and time for me ZZZ",
    "And time yet for a hundred indecisions XXX",
    "And for a hundred visions and revisions XXX",
    "Before the taking of a toast and tea ZZZ.",
    "In the room the women come and go XXX",
    "Talking of Michelangelo. XXX",
    "And indeed there will be time XXX",
    'To wonder, "Do I dare?" and, "Do I dare?" ZZZ',
    "Time to turn back and descend the stair, ZZZ",
    "With a bald spot in the middle of my hair — XXX",
    '(They will say: "How his hair is growing thin!") XXX',
    "My morning coat, my collar mounting firmly to the chin, ZZZ",
    "My necktie rich and modest, but asserted by a simple pin — XXX",
    '(They will say: "But how his arms and legs are thin!") ZZZ',
    "Do I dare XXX",
    "Disturb the universe? XXX",
    "In a minute there is time ZZZ",
    "For decisions and revisions which a minute will reverse. XXX",
]


targets = [0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1,
           1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1]
```

## Naive Bayes from scratch

> The focus of this section is to get the maths right -- we don't care much about the efficiency or elegancy of the code.

The first thing we should think about is the data we're feeding into the algorithm. The data at hand are sentences, so we need to parse them into words or tokens.

### The tokenizer

We'll take a very simple approach here -- given a sentence, our function returns a dictionary with words in the sentence as keys and counts of the words as values[^tokenizer]. We also define a helper function that adds counts to our result dictionaries.

[^tokenizer]: Here we only keep tokens that are at least 2 characters long to stay consistent with the `sklearn.feature_extraction.text.CountVectorizer` method.

```python
import re

def add_count(x, output_dict, count=1):
    if x not in output_dict:
        output_dict[x] = count
    else:
        output_dict[x] += count

def tokenize(sentence):
    tokens = {}
    # Only keep tokens that are at least 2 characters long
    for token in re.findall(r"\b\w\w+\b", sentence):
        token = token.lower()
        add_count(token, tokens)
    return tokens
```

### Naive Bayes core

Let's take a minute to think about the parts needed to calculate the final probability for each class:

1. For each token, find its count in each of the classes.
2. Calculate and store $P(C)$, the frequencies of the classes.
3. Calculate the (log) probability of each token for each class using Bayes' Theorem. Laplace smoothing should be applied.
4. Predict the class of a given sentence by returning the class with maximum probability.

Let's get started! First we take down everything that needs to be stored:

```python
class MyNaiveBayes:
    def __init__(self, laplace_smoothing_param=1):
        self.k = laplace_smoothing_param
        self.classes_count = {}
        self.class_tokens = {}
        self.all_tokens = set()
```

We want to keep the class flexible and not limit the input to two classes. `classes_count` would have different classes as keys and their counts as values. `class_tokens` would be a nested dictionary -- each element is a dictionary of token counts with the key as the class label.

Training the model is pretty easy as we only need to go through the corpus once. We first go through each sentence and count the words. Then we collect all the observed words into `all_tokens` and add the Laplace pseudo-counts to `class_tokens`.

```python
def train(self, corpus, labels):
    for label, s in zip(labels, corpus):
        add_count(label, self.classes_count)
        if label not in self.class_tokens:
            self.class_tokens[label] = {}

        tokens = tokenize(s)
        self.all_tokens.update(tokens.keys())
        for token, token_count in tokens.items():
            add_count(token, self.class_tokens[label], token_count)

    for label in self.class_tokens.keys():
        for token in self.all_tokens:
            if token not in self.class_tokens[label]:
                self.class_tokens[label][token] = 0
            self.class_tokens[label][token] += self.k  # laplace
```

Calculating the probability for each class is trivial:

```python
def calc_class_proba(self, classID):
    return self.classes_count[classID] / np.sum(
        [x for x in self.classes_count.values()]
    )
```

Finding the probability for each token given the class is only slightly harder -- we run into cases where the word count is zero, so the Laplace estimator is plugged in:

```python
def calc_token_proba(self, token, classID):
    all_tokens_in_class = np.sum([x for x in self.class_tokens[classID].values()])
    return self.class_tokens[classID][token] / all_tokens_in_class
```

To get the probability for each class, we just apply Bayes' Theorem and scale the sum of the probabilities to one.

```python
def predict_proba(self, sentence):
    tokens = tokenize(sentence)
    res = {}
    for label in self.classes_count.keys():
        # start from the class probability
        log_proba = np.log2(self.calc_class_proba(label))
        for token, token_count in tokens.items():
            # handle unseen words
            if token in self.class_tokens[label]:
                token_proba = self.calc_token_proba(token, label)
                log_proba += token_count * np.log2(token_proba)

        res[label] = np.exp2(log_proba)

    denom = np.sum([x for x in res.values()])
    res = {label: proba / denom for label, proba in res.items()}
    return res
```

And that's pretty much the core part of our model! We're just missing the functions to predict the final class, which is given below in the full class with all the other methods.

```python
import re

import numpy as np

def add_count(x, output_dict, count=1):
    if x not in output_dict:
        output_dict[x] = count
    else:
        output_dict[x] += count


def tokenize(sentence):
    tokens = {}
    # Only keep tokens that are at least 2 characters long
    for token in re.findall(r"\b\w\w+\b", sentence):
        token = token.lower()
        add_count(token, tokens)
    return tokens


class MyNaiveBayes:
    def __init__(self, laplace_smoothing_param=1):
        self.k = laplace_smoothing_param
        self.classes_count = {}
        self.class_tokens = {}
        self.all_tokens = set()

    def train(self, corpus, labels):
        for label, s in zip(labels, corpus):
            add_count(label, self.classes_count)
            if label not in self.class_tokens:
                self.class_tokens[label] = {}

            tokens = tokenize(s)
            self.all_tokens.update(tokens.keys())
            for token, token_count in tokens.items():
                add_count(token, self.class_tokens[label], token_count)

        for label in self.class_tokens.keys():
            for token in self.all_tokens:
                if token not in self.class_tokens[label]:
                    self.class_tokens[label][token] = 0
                self.class_tokens[label][token] += self.k  # laplace

    def calc_class_proba(self, classID):
        return self.classes_count[classID] / np.sum(
            [x for x in self.classes_count.values()]
        )

    def calc_token_proba(self, token, classID):
        all_tokens_in_class = np.sum([x for x in self.class_tokens[classID].values()])

        # handle unseen words
        if token not in self.class_tokens[classID]:
            return 1
        return self.class_tokens[classID][token] / all_tokens_in_class

    def predict_proba(self, sentence):
        tokens = tokenize(sentence)
        res = {}
        for label in self.classes_count.keys():
            # start from the class probability
            log_proba = np.log2(self.calc_class_proba(label))
            for token, token_count in tokens.items():
                # handle unseen words
                if token in self.class_tokens[label]:
                    token_proba = self.calc_token_proba(token, label)
                    log_proba += token_count * np.log2(token_proba)

            res[label] = np.exp2(log_proba)

        denom = np.sum([x for x in res.values()])
        res = {label: proba / denom for label, proba in res.items()}
        return res

    def predict(self, sentence):
        proba = self.predict_proba(sentence)
        max_proba = max([x for x in proba.values()])
        res = [x for x in proba.keys() if proba[x] == max_proba]
        return res[0]
```

### Checking the classifier

Now we check our model with the corpus given at the beginning. We have a test set given below:

```python
new_data = [
    "For I have known them all already, known them all: ZZZ",
    "Have known the evenings, mornings, afternoons, ZZZ",
    "I have measured out my life with coffee spoons; XXX",
    "I know the voices dying with a dying fall XXX",
    "Beneath the music from a farther room. ZZZ",
    "So how should I presume? ZZZ",
]
actual = [0, 0, 1, 1, 0, 0]

cls = MyNaiveBayes()
cls.train(corpus, targets)
[cls.predict_proba(x) for x in new_data]
# [{0: 0.9721513447351029, 1: 0.027848655264897063},
#  {0: 0.9026159391741999, 1: 0.09738406082580003},
#  {0: 0.02876761464981322, 1: 0.9712323853501869},
#  {0: 0.021732005477143462, 1: 0.9782679945228566},
#  {0: 0.8132271179857358, 1: 0.18677288201426423},
#  {0: 0.9251393970016891, 1: 0.07486060299831092}]
```

Well, we got all of the test cases right! Next we'll compare our results to the output of `scikit-learn` to see if the probabilities are correct.

## Using scikit-learn

The `CounterVectorizer` class in `scikit-learn` can be used to convert each string in the corpus into a vector of word counts. The output is a compressed sparse row matrix (`scipy.sparse.csr_matrix`)[^count-vectorizer-output].

```python
from sklearn.feature_extraction.text import CountVectorizer
count_vect = CountVectorizer()
X_train_counts = count_vect.fit_transform(corpus.X)
```

[^count-vectorizer-output]:
    We can check the data stored in the sparse matrix by converting it to a pandas data frame:

    ```python
    names = count_vect.get_feature_names()
    arr = X_train_counts.toarray()
    counts_df = pd.DataFrame(arr, columns=names)
    ```

Once we have the vectorized input, we can create a naive Bayes classifier by fitting it to the corpus. Then we may invoke the `predict` method to classify new instances.

```python
from sklearn.naive_bayes import MultinomialNB

clf = MultinomialNB().fit(X_train_counts, targets)
clf.predict_proba(X_new_counts)
# array([[0.97215134, 0.02784866],
#        [0.90261594, 0.09738406],
#        [0.02876761, 0.97123239],
#        [0.02173201, 0.97826799],
#        [0.81322712, 0.18677288],
#        [0.9251394 , 0.0748606 ]])
```

The output values are almost identical to that of our hand-written model above.

## Introducing NLTK

[NLTK](https://www.nltk.org/) (the Natural Language ToolKit) provides a suite of text processing libraries for classification, tokenization, stemming, tagging, etc. It also provides interfaces to over 50 corpora and lexical resources such as WordNet. Here we'll briefly introduce some of the most common methods. For detailed tutorials, see [this book by Steven Bird et al.](http://www.nltk.org/book/).

### Tokenization

NLTK has tokenizers for splitting text into sentences (based on capitalization and punctuation) and into words. They should be much more robust than the one we used above.

```python
import nltk
from nltk.corpus import wordnet
from nltk.stem import PorterStemmer, WordNetLemmatizer
from nltk.tokenize import sent_tokenize, word_tokenize

poem = ". ".join(corpus)
sentences = sent_tokenize(poem)
print(sentences[:5])
# ['And indeed there will be time ZZZ.',
#  'For the yellow smoke that slides along the street XXX.',
#  'Rubbing its back upon the window-panes ZZZ.',
#  'There will be time, there will be time ZZZ.',
#  'To prepare a face to meet the faces that you meet XXX.']

words = word_tokenize(poem)
print(words[:5])
# ['And', 'indeed', 'there', 'will', 'be']
```

### Stemming and lemmatization

As defined by the Stanford NLP group, the goal of both stemming and lemmatization is to reduce inflectional forms and sometimes derivationally related forms of a word to a common base form. `Stemming` usually refers to a crude heuristic process that chops off the ends of words in the hope of achieving this goal correctly most of the time, and often includes the removal of derivational affixes. `Lemmatization` usually refers to doing things properly with the use of a vocabulary and morphological analysis of words, normally aiming to remove inflectional endings only and to return the base or dictionary form of a word.

Let's see the stemmed word tokens first. All characters are converted to lowercase, and some errors were introduced, e.g. "quickly" to "quickli".

```python
text = "Run runs running ran. Tall tallest. Quick quickest quickly."
words = word_tokenize(text)

stemmer = PorterStemmer()
stemmed = [(tok, stemmer.stem(tok)) for tok in words]
print(stemmed)
# [('Run', 'run'),
#  ('runs', 'run'),
#  ('running', 'run'),
#  ('ran', 'ran'),
#  ('.', '.'),
#  ('Tall', 'tall'),
#  ('tallest', 'tallest'),
#  ('.', '.'),
#  ('Quick', 'quick'),
#  ('quickest', 'quickest'),
#  ('quickly', 'quickli'),
#  ('.', '.')]
```

We use WordNet Synset to demonstrate lemmatization. WordNet is a lexical database designed for NLP in English, and Synset is an interface to look up words in WordNet. It takes slightly longer to run, but the results are more accurate.

```python
lem = WordNetLemmatizer()
lemmatized = [(tok, lem.lemmatize(tok)) for tok in words]
print(lemmatized)
# [('Run', 'Run'),
#  ('runs', 'run'),
#  ('running', 'running'),
#  ('ran', 'ran'),
#  ('.', '.'),
#  ('Tall', 'Tall'),
#  ('tallest', 'tallest'),
#  ('.', '.'),
#  ('Quick', 'Quick'),
#  ('quickest', 'quickest'),
#  ('quickly', 'quickly'),
#  ('.', '.')]
```

### Part of speech tagging

Another common use case of NLTK is to tag the part of speech[^part-of-speech] given a list of tokens. This is useful for identifying entities and relationships between entities in the text.

[^part-of-speech]: See `nltk.help.upenn_tagset()` for all of the tags.

```python
tokens = word_tokenize("The very first run was unsuccessful. ")
tagged = nltk.pos_tag(tokens)
for pair in tagged:
    print(pair)
    nltk.help.upenn_tagset(pair[1])

# ('The', 'DT')
# DT: determiner
#     all an another any both del each either every half la many much nary
#     neither no some such that the them these this those
# ('very', 'RB')
# RB: adverb
#     occasionally unabatingly maddeningly adventurously professedly
#     stirringly prominently technologically magisterially predominately
#     swiftly fiscally pitilessly ...
# ('first', 'JJ')
# JJ: adjective or numeral, ordinal
#     third ill-mannered pre-war regrettable oiled calamitous first separable
#     ectoplasmic battery-powered participatory fourth still-to-be-named
#     multilingual multi-disciplinary ...
# ('run', 'NN')
# NN: noun, common, singular or mass
#     common-carrier cabbage knuckle-duster Casino afghan shed thermostat
#     investment slide humour falloff slick wind hyena override subhumanity
#     machinist ...
# ('was', 'VBD')
# VBD: verb, past tense
#     dipped pleaded swiped regummed soaked tidied convened halted registered
#     cushioned exacted snubbed strode aimed adopted belied figgered
#     speculated wore appreciated contemplated ...
# ('unsuccessful', 'JJ')
# JJ: adjective or numeral, ordinal
#     third ill-mannered pre-war regrettable oiled calamitous first separable
#     ectoplasmic battery-powered participatory fourth still-to-be-named
#     multilingual multi-disciplinary ...
# ('.', '.')
# .: sentence terminator
#     . ! ?
```
