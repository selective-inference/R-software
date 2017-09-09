Rcpp: 
	- rm -f selectiveInference/src/RcppExports.cpp
	- rm -f selectiveInference/R/RcppExports.R
	Rscript -e "library(Rcpp); Rcpp::compileAttributes('selectiveInference')"

install: Rcpp
	R CMD INSTALL selectiveInference

build: 
	R CMD build selectiveInference

check: Rcpp build
	R CMD build selectiveInference
	R CMD check selectiveInference_1.2.2.tar.gz # fix this to be a script variable