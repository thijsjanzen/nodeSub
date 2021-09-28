# nodeSub <img src="pics/nodesub_sticker.png" align="right" width="200" />

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/nodeSub)](https://cran.r-project.org/package=nodeSub)
[![R-CMD-check](https://github.com/thijsjanzen/nodeSub/workflows/R-CMD-check/badge.svg)](https://github.com/thijsjanzen/nodeSub/actions)

Branch |[![Codecov logo](pics/Codecov.png)]
---|---
master |[![codecov.io](https://codecov.io/gh/thijsjanzen/nodeSub/branch/master/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/nodeSub/)
develop|[![codecov.io](https://codecov.io/gh/thijsjanzen/nodeSub/branch/develop/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/nodeSub/)
richel |[![codecov.io](https://codecov.io/gh/thijsjanzen/nodeSub/branch/richel/graph/badge.svg)](https://codecov.io/gh/thijsjanzen/nodeSub/)

NodeSub is an R package that can be used to generate alignments under the node substitution model.

# Installation

NodeSub also provides functions to perform BEAST2 analyses, and if you want access to those, you will need to rely on the functionality of the babette suite. Unfortunately, the babette suite is currently not available on CRAN, but can be installed from rOpenSci:

```
remotes::install_github("ropensci/beautier")
remotes::install_github("ropensci/tracerer")
remotes::install_github("ropensci/beastier")
remotes::install_github("ropensci/mauricer")
remotes::install_github("ropensci/babette")
remotes::install_github("ropensci/mcbette")
```
Next, we need to install BEAST2 and the BEAST2 NS package (for marginal likelihood calculation) - and make sure R can find it. This we can do as follows:

```
remotes::install_github("richelbilderbeek/beastierinstall") 
beastierinstall::install_beast2() 
beastier::is_beast2_installed()
remotes::install_github("richelbilderbeek/mauricerinstall") 
mauricerinstall::install_beast2_pkg("NS")
mauricer::is_beast2_ns_pkg_installed()
```

# Usage
Given a tree, we can generate a NodeSubstitution alignment as follows:
```
input_tree <- TreeSim::sim.bd.taxa(n = 100,
                                   numbsim = 1,
                                   lambda = 1,
                                   mu = 0.1,
                                   complete = TRUE)[[1]]
                                     
target_alignment <- sim_unlinked(phy = input_tree,
                                 rate1 = 1e-3,
                                 rate2 = 1e-3,
                                 l = 10000,
                                 node_time = 0.3)   
```
We can then proceed to create a Twin alignment (e.g. an alignment with the exact same number of accumulated substitutions, given the same tree, but using a 'normal' substitution model)
```
comp_alignment <- create_equal_alignment(input_tree = geiger::drop.extinct(input_tree),  # can only work on trees without extinct branches
                                         sub_rate = 1e-3,
                                         alignment_result = target_alignment)
```
Now we have two alignments, one generated using the Node Substitution model, and one using standard strict clock model. For both we can perform phylogenetic inference to get the resulting tree. The aim is to compare the posterior distribution of trees with the true tree that we started with 'input_tree', and estimate the error invoked by the node substitution model.

```
node_posterior <- infer_phylogeny(target_alignment$alignment,
                                  "node_posterior",
                                  clock_prior = beautier::create_strict_clock_model(clock_rate_param = beautier::create_clock_rate_param(value = 1e-3)),
                                  burnin = 0.1,
                                  working_dir = get_wd(),
                                  sub_rate = 1e-3)
                                  
reference_posterior <- infer_phylogeny(comp_alignment$alignment,
                                       "reference_posterior",
                                       burnin = 0.1,
                                       clock_prior = beautier::create_strict_clock_model(clock_rate_param = beautier::create_clock_rate_param(value = comp_alignment$adjusted_rate)),
                                       working_dir = getwd(),
                                       sub_rate = 1e-3)                               
```
Having two posterior distributions of trees, we can compare them based on summary statistics.
```
node_stats <- calc_sum_stats(node_posterior$all_trees,
                             input_tree)
ref_stats  <- calc_sum_stats(reference_posterior$all_trees,
                            input_tree)
```                                      
