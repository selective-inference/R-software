Rcpp: 
	rm selectiveInference/src/RcppExports.cpp
	rm selectiveInference/R/RcppExports.R
	Rscript -e "library(Rcpp); Rcpp::compileAttributes('selectiveInference')"