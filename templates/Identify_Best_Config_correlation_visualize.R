# add libraries
library(dplyr) # includes >%> and group_by
library(purrr) # includes map (https://rebeccabarter.com/blog/2019-08-19_purrr, https://github.com/rstudio/cheatsheets/blob/main/purrr.pdf), map 2, and (i)walk
library(yaml) # includes read_yaml and write_yaml

# set working directory
setwd("C:/Users/slfri/ownCloud/7T_WISC_MVPA/derivatives/correlation/grOWL/visualize")
# set path to where data are stored on CHTC
data_path <- '/staging/groups/psych_rogers_group/saskia/7T/data/wholebrain/averaged'

# read .yaml used for setting up tune directory structure
stub_tune <- list(read_yaml("tune.yaml"))
# if you want to do several at once, use this setup
# stub_tune <- list(
#     SIFT = read_yaml("./SIFT/correlation/performance/tune/stub_sift.yaml"),
#     CAFFE7 = read_yaml("./CAFFE7/correlation/performance/tune/stub_caffe7.yaml"),
#     CAFFE8 = read_yaml("./CAFFE8/correlation/performance/tune/stub_caffe8.yaml"),
#     NEXT6D = read_yaml("./NEXT6D/correlation/performance/tune/stub_next6D.yaml")
# )

# load .csv generated from condensing directory tree
tune_df <- list(read.csv("tune_performance.csv"))
# if you want to do several at once, use this setup
# tune_df <- list(
#     SIFT = read_csv("./SIFT_growl_correlation_performance_tune.csv"),
#     CAFFE7 = read_csv("./CAFFE7_growl_correlation_performance_tune.csv"),
#     CAFFE8 = read_csv("./CAFFE8_growl_correlation_performance_tune.csv"),
#     NEXT6D = read_csv("./NEXT6D_growl_correlation_performance_tune.csv")
# )

# get the model's loss on the test set (err1), the model's loss on the training
# set (err2) and the number of nonzero voxels. 
# To achieve this, first group the rows into groups with the same data, metadata, 
# final holdout set, and combinations of lambda and lambda1. Then find the mean
# test-set loss (err1), training set loss (err2) and number of nonzero voxels for
# each group. Get rid of columns that we don't need. Sort the rows into 
# groups with the same data and final holdout, Then find, for each group, 
# the row with the lowest test-set loss (i.e. the combination of lambda and lambda1
# that works best for that participant based on those 9 folds of training data).
# Finally, get rid of groups.
# Notes:
# - %>% means "and then" (and is an alternative to nesting functions)
# - ~{.x} is short for function(x){}
# - we use map so that we can store data from multiple .csvs in a single list
#   and the function can be applied to the list (rather than to a data frame)

best_df <- tune_df %>% map(., ~{
  .x %>%
    group_by(metadata, data, finalholdout, lambda, lambda1)  %>%
    summarize(
      err1 = mean(err1),
      err2 = mean(err2),
      nzvox = mean(nzvox)
    ) %>%
    group_by(data, finalholdout)  %>%
    slice_min(err1, n = 1) %>%
    ungroup()
})

# adjust the contents of stub_tune, using the contents of best_df, to make a list of all the info
# needed to make the .yaml for the final phase
stub_final <- map2(stub_tune, best_df, ~{
  .x$metadata <- file.path(data_path, .y$metadata)
  .x$data <- file.path(data_path, .y$data)
  .x$finalholdout <- 0L # forces 0 to be a integer
  .x$cvholdout <- as.integer(.y$finalholdout[1])
  .x$lambda <- .y$lambda
  .x$lambda1 <- .y$lambda1
  .x$SearchWithHyperband <- FALSE
  .x$HYPERBAND <- NULL
  .x$BRACKETS <- NULL
  .x$EXPAND <- list(c('data','metadata','lambda','lambda1'))
  return(.x)
})

# adjust that list to make a list of all the info needed to make the .yaml for the perm phase
# N.B. the line beginning x$RandomSeed describes how many permutations will be assigned to each job. 
# split(1:A, rep(1:B, each = C)) splits a total of A permutations (A/10 per fold) into B jobs and each
# job will do C permutations.
# N.B. the total number of jobs is B x 10 (folds) x n participants. B should be chosen such that
# the total number of jobs should be under 10,000, otherwise condor_q will not be able to queue them.
# N.B. There should be at least 1000 permutations per participant unless you are planning on computing
# a permuted-group-average.
stub_perm <- stub_final %>%
  map(~{
    .x$EXPAND <- list(c('data','metadata','lambda','lambda1'),'RandomSeed')
    .x$URLS <- c('data','metadata','PermutationIndex')
    .x$PermutationIndex <- file.path(dirname(.x$metadata[1]), "PERMUTATION_STRUCT.mat")
    .x$PermutationTest <- TRUE
    .x$PermutationMethod <- "manual"
    .x$perm_varname <- "PERMUTATION_STRUCT"
    .x$RandomSeed <- unname(split(1:100, rep(1:10, each = 10)))
    return(.x)
  })


# write final .yaml.
# N.B. walk is similar to map, but is used for writing outputs
stub_final %>%
  walk(~{
    write_yaml(.x, file = sprintf('final.yaml'))
  })

# write perm .yaml
stub_perm %>%
  walk(~{
    write_yaml(.x, file = sprintf('perm.yaml'))
  })

# if you want to do several at once, use this setup
# stub_final %>%
#   iwalk(~{
#     write_yaml(.x, file = sprintf('stub_final_%s_correlation_performance.yaml', .y))
#   })
# 
# stub_perm %>%
#   iwalk(~{
#     write_yaml(.x, file = sprintf('stub_perm_%s_correlation_performance.yaml', .y))
#   })