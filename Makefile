Rcpp: 
	- rm -f selectiveInference/src/RcppExports.cpp
	- rm -f selectiveInference/R/RcppExports.R
	Rscript -e "library(Rcpp); Rcpp::compileAttributes('selectiveInference')"