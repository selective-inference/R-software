Rcpp: 
	- rm -f selectiveInference/src/RcppExports.cpp
	- rm -f selectiveInference/R/RcppExports.R
	- Rscript -e "library(Rcpp); Rcpp::compileAttributes('selectiveInference')"

install: Rcpp src
	R CMD INSTALL selectiveInference

build: src 
	R CMD build selectiveInference

src:
	cp C-software/src/* selectiveInference/src

check: Rcpp build 
	R CMD build selectiveInference
	R CMD check selectiveInference_1.2.5.tar.gz --as-cran # fix this to be a script variable