# Sparse Overlapping Sets (SOS) Lasso

## Unfiled note to self!
The `subject_id_fmt` field, and ultimately `sscanf` that uses it, can be
used to extract a wide variety of substrings. Notably, if it is
extracting only a number, it will return the value as a numeric data
type, but in other cases it is perfectly happy to return a string.

Consider the following examples (where the second argument is the value
of `subject_id_fmt`. In the first, we use `sscanf` to extract the digits
`101` from the filename and return them as the number `101`. In the
second, rather than matching digits, we match an unbroken string of
characters that are each in the range 0-9, and return a string. Finally,
by including 'BM' within the set of characters to match, we can extract
the full subject code 'BM101', and trim the (in this case) unwanted
information. (Note, the goal is to isolate the portion of the filename
that corresponds to the `subject` field within the `metadata` structure.
This is how metadata and data are related to one another at present).

```
sscanf('BM101_avg.mat','BM%d_avg.mat')
ans =
    101
```

```
sscanf('BM101_avg.mat','BM%[0-9]_avg.mat')
ans =
'101'
```

```
sscanf('BM101_avg.mat','%[BM0-9]_avg.mat')
ans =
'BM101'
```

Note that this last option is not very precise. It would match any
permutation of the numbers 0-9 and the letters B and M. If precision is
important, something like this would match exactly BM[numbers].

```
sscanf('BM101_avg.mat','%[B]%[M]%[0-9]_avg.mat')
```

## Overview and background

### Standard regression
Regression is a way of learning a model that will make predictions
about some variable of interest based on the states of several
features. The resulting model is a vector of weights, one per feature.
A limitation of regression is that, in order to find a well
defined model, there must be more observations in the training set
than there are features.

In general, regression is minimizing the difference between some
observed "ground truth", `y`, and the weighted sum of some features
that that potentially covary with `y` in some meaningful way, `X`.
The weight vector itself is noted as `b`. So, roughly speaking, the
goal is to minimize:

    f(y,X*b)

Where `f()` is a loss function, such as "squared error" if `y` is
continuous or "logistic loss" if `y` is binary. `X` is a matrix,
where columns are features and rows are observations, and `b` is the
vector of weights that is being learned in order to minimize `f()`.
`X*b` is the inner product of `X` and `b`.

### When standard regression falls short
In many cases, there are many more features than observations, and
the objective is to identify the features that are actually useful in
the context of a given classification problem. For example, cognitive
neuroscientists are interested in the patterns of neural activity
that are associated with different mental states. Given an fMRI
dataset, which may have tens or hundreds of thousands of features
(i.e., voxels) associated with each observation (i.e., TR, trial,
condition, etc.). Critically, there will always be more features than
observations unless the dataset is somehow reduced.

### Lasso (least absolute shrinkage and selection operator)
There are many ways to go about reducing the number of features.
Here, we consider a solution called Lasso [1]. Lasso involves
modifying the regression problem. In addition simply minimizing `f()`,
an additional component is added to the optimization:

    f(y,X*b) + lambda*h(b)

Now, the goal is to *jointly* minimize `f()` and `h()`, where `h()`
is a function ONLY of the weights. That is, it does not have to do
with prediction, but the weights themselves. In the case of Lasso,
`h()` is the sum of the absolute values of the weights. This means
that `h()` is minimized when all weights are 0. Of course, this will
not make for good predictions. The goal is the find a weight vector
`b` that obtains good prediction accuracy (i.e., a small value of
`f()`) while utilizing as few features as possible. How much `h()`
matters relative to `f()` can be adjusted using the free parameter
`lambda`.

### Limits of Lasso
While on one hand Lasso solves the problem of having more features
than observations by seeking a "sparse" solution, it obtains sparsity
by selecting a set of features that each contributes as much unique
information as possible. That is, Lasso will provide a solution where
the selected features are as non-redundant as possible. In other
words, Lasso may be **highly selective** in that the voxels it
identifies tend to carry meaningful signal, but have **low
sensitivity** in that it will tend to under-report the number of features
that are truly relevant to the problem.

It is also important to keep in mind the information available to
Lasso as it seeks an optimal weight vector `b`. All it knows is that
there are a set of features, and the states of those features might
be used to make accurate predictions that can minimize `f()`. If
those features have spatial relations, or if groups of these features
should be considered together, Lasso knows nothing about that. The
only information that contributes to the weight Lasso assigns to
a particular feature are the values returned from `f()` and `h()`.

Finally, Lasso is applied to only one dataset at a time. Often,
several datasets may bear on a common question, and one might expect
that similar features ought to be discovered in each. Lasso is
incapable of taking advantage of the structure across datasets. In
terms of fMRI data, Lasso is an inherently single-subject analysis
approach.

### SOS Lasso
SOS Lasso [2] attempts to address these limitations by allowing features
to be organized into sets. Instead of penalizing all weights equally,
as is done by Lasso, SOS Lasso will penalize less within sets than
across sets. It involves adding yet another component to the
optimization:

    f(y,X*b) + lambda*h(b) + gamma*g(s,b)

Where `s` contains the set label for each feature. Whereas the Lasso
penalty `h()` evaluates to the sum of the absolute value of the
weights in `b`, `g()` effectively loops over groups, takes the square
root of the sum of the weights within each group, and then sums over
those.

The formulation above was simply illustrative, and actually differs from
the implementation in an important respect. Rather than having a totally
independent free parameter on each component, we have formulated the
penalty so that `lambda` scales the sum of `h()` and `g()`, and a second
parameter `alpha` titrates between the importance of `h()` vs. `g()`:

    f(y,X*b) + lambda*( (1-alpha*h(b)) + alpha*g(s,b) )

This has the advantage of bounding `alpha` and making it a little more
interpretable: if `alpha=1`, then the value of `h()` is set to zero and
the only thing that matters is `g()`, the "sets penalty", while if
`alpha=0`, the optimization reduces to Lasso (because `g()` will be set
to zero).

#### Prefering features from the same set
Consider a simple case where there are 2 groups, and 2 alternative
`b` vectors. Each has two non-zero weights, where the non-zero values
are 1. In `b1`, the ones belong to different groups, and in `b2` they
belong to the same group. The following (purely illustrative) Matlab
code demonstrates how `g()` evalutes to a lower value in the  the
"more structured" case.

```matlab
s  = [1, 1, 2, 2]
b1 = [0, 1, 0, 1]
b2 = [1, 1, 0, 0]

gval1 = zeros(1,2);
gval2 = zeros(1,2);

for i = 1:2
  gval1(i) = sqrt(sum(b1(s==i)));
  gval2(i) = sqrt(sum(b2(s==i)));
end

sum(gval1)
sum(gval2)
```

#### More sensitive
Since features from the same set can be utilized more cheaply than
an arbitrary set of features, this has the potential to increase
sensitivity if the true signal really does belong to the defined
sets.

#### Multi-task learning
SOS Lasso also allows for multi-task learning. In short, this means
that instead of arriving at a completely independent solution for
each dataset as would be the case with standard Lasso, it is possible
to fit several models simultaneously whose solutions are dependent on
one another. This is achieved by using the same group labels across
datasets. In reality when `g()` is evaluated, it sums over the
weights within a group **across all datasets**. So, if a feature is
selected from a particular set in one solution, it will be "easier"
or "less costly" to select a unit from that same group in a different
dataset.

### Important points of note

#### OVERLAPPING Sets
It is possible for a feature to belong to more than one set.
Practically, this means that you can make an assumption like "If a
feature is informative, that increases the likelihood that a feature
that corresponds to a neighboring point in space will also be
informative" without making strict commitments about how the space
should be carved up. This means you do not need to approach the
problem with a strong hypothesis about how the features ought to be
grouped.

#### SPARSE sets
SOS Lasso does not select sets in their entirety. Assigning a weight
to one feature does NOT necessitate assigning a weight to all other
members of the set. SOS Lasso simply prefers to select features from
a small number of sets (see the illustrative example above). Sparsity
is achieved even within sets.

## References
[1] Tibshirani (1996). "Regression shrinkage and selection via
the lasso". Journal of the Royal Statistical Society, Series B 58
(1): 267–288

[2] Rao, Cox, Nowak, and Rogers (2013). "Sparse Overlapping Sets Lasso
for multitask learning and its application to fMRI analysis".
Advances in Neural Information Processing Systems 26, 2202--2210
