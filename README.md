# R-software
R software for selective inference
Authors: Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid
Maintainer: Rob Tibshirani <tibs@stanford.edu>

New tools for inference after selection, for use with forward stepwise regression, least angle regression, the lasso, and the many means problem.

Code is in the directory selectiveInference/R.
* funs.common.R: Basic functions used by many other functions, such as standardization.
* funs.fixed.R: Inference for LASSO at a fixed, deterministic value of lambda.
* funs.fs.R: Inference for forward stepwise.
* funs.groupfs.R: Inference for forward stepwise with groups of variables, e.g. factors.
* funs.inf.R: Common functions for inference with fixed, fs, lar, and manymeans (but not group).
* funs.lar.R: Inference for least angle regression.
* funs.max.R: Some numerical approximations. Deprecated?

