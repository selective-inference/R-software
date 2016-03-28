# R-software
R software for [selective inference](http://cran.r-project.org/web/packages/selectiveInference/).  
Authors: Ryan Tibshirani, Rob Tibshirani, Jonathan Taylor, Joshua Loftus, Stephen Reid  
Maintainer: Rob Tibshirani <tibs@stanford.edu>

New tools for inference after selection, for use with forward stepwise regression, least angle regression, the lasso, and the many means problem. The package is available on [CRAN](http://cran.r-project.org/web/packages/selectiveInference/). See [this paper](http://www.pnas.org/content/112/25/7629.full) for a high level introduction to selective inference.

Code is in the directory selectiveInference/R.
* funs.common.R: Basic functions used by many other functions, such as standardization.
* funs.fixed.R: Inference for LASSO at a fixed, deterministic value of lambda.
* funs.fs.R: Inference for forward stepwise.
* funs.groupfs.R: Inference for forward stepwise with groups of variables, e.g. factors.
* funs.inf.R: Common functions for inference with fixed, fs, lar, and manymeans (but not group).
* funs.lar.R: Inference for least angle regression.
* funs.max.R: Some numerical approximations. Deprecated?
