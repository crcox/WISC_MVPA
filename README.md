[![DOI](https://zenodo.org/badge/19497385.svg)](https://zenodo.org/badge/latestdoi/19497385)

# Whole-brain Imaging with Sparse Correlations (WISC) Workflow

The WISC workflow implements whole-brain classification and representational analysis methods for cognitive neuroscience. It highlights two techniques developed at the University of Wisconsin-Madison and the Wisconsin Institutes for Discovery:

1. Sparse Overlapping Sets (SOS) Lasso [1]
2. Network Representational Similarity Analsysis [2]

The tool is intended for use with Hight Throughput/Performance Computing clusters: parameters can be specified using a JSON formatted text file in the working directory, and results are written to disk as a structure. Tools have been written to aggregate over `results.mat` files returned by several parallel jobs.

## Usage
Adapt the `Makefile` to point to your (or your institution's) Matlab Compiler `mcc`. This will build a portable executable that can be run anywhere there is a compatible Matlab Runtime Environment without duplicating licenses. Your computing environment may have specific documentation on this process. I build the code against **Matlab 2015b**, and avoid using data types and features introduced in recent years. However, Matlab's feature set changes all the time time making forward or backward compatability difficult to predict.

Once built, running the executable in the same directory as a JSON file called `params.json` will instruct the program. JSON fields should correspond to the parameters listed at the beginning on `src/WISC_MVPA.m`.

When completed, the program will return a file called `results.mat`. This will contain a structured array that will differ depending on whether you are running **SOS LASSO** or **Network RSA** (aka GrOWL).

## Data
Functional imaging data are assumed to be represented over two Matlab `.mat` files. The first contains the functional imagine data itself for a single subject. This file will contain one or more matrices (if there are multiple versions of the data for each subject). Which matrix is loaded is controlled with the `data_varname` parameter, which is a required parameter (even if only one matrix exists in the file). Rows of this matrix correspond to training examples (trials, stimuli, etc.) and columns correspond to features (voxels, electrodes, etc.).

The second file contains metadata. Metadata corresponds to targets, filters, coordinates, crossvalidation labels, and more. See the `demos/` folder for examples defining data and metadata files and for running basic analyses on simulated data.

Paths to data (and variables to read from these `.mat` files) is all specified within the params.json file. 

# Sparse Overlapping Sets (SOS) Lasso

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
Here, we consider a solution called Lasso [3]. Lasso involves
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
SOS Lasso [1] attempts to address these limitations by allowing features
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

# Network RSA
Network RSA leverages the concepts of structured sparsity and multitask learning in a different way. Rather than multiple subjects related to a 1-dimensional target structure (a classification problem), Network RSA relates the data from one subject to a multidimensional target structure.

Many of the concepts covered above are relevant to Network RSA as well. For a full discussion of the technique please see our publication [2].

## References
[1] Rao, Cox, Nowak, and Rogers (2013). "Sparse Overlapping Sets Lasso
for multitask learning and its application to fMRI analysis".
Advances in Neural Information Processing Systems 26, pp 2202--2210.

[2] Oswal, Cox, Lambon Ralph, Rogers, Nowak (2016). "Representational similarity learning with application to brain networks". ICML'16 Proceedings of the 33rd International Conference on Machine Learning, 48, pp 1041--1049.

[3] Tibshirani (1996). "Regression shrinkage and selection via the lasso". Journal of the Royal Statistical Society, Series B 58(1), pp 267-–288.
