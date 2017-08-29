#include <Rcpp.h>      // need to include the main Rcpp header file 
#include <matrixcomps.h>    // where update1, downdate1 are defined

// [[Rcpp::export]]
Rcpp::List update1_(Rcpp::NumericMatrix Q2,
		    Rcpp::NumericVector w,
		    int m,
		    int k) {

  update1(Q2.begin(),
	  w.begin(),
	  m,
	  k);

  return(Rcpp::List::create(Rcpp::Named("Q2") = Q2,
			    Rcpp::Named("w") = w));
}

// [[Rcpp::export]]
Rcpp::List downdate1_(Rcpp::NumericMatrix Q1,
		      Rcpp::NumericMatrix R,
		      int j0,
		      int m,
		      int n) {

  downdate1(Q1.begin(),
	    R.begin(),
	    j0,
	    m,
	    n);

  return(Rcpp::List::create(Rcpp::Named("Q1") = Q1,
			    Rcpp::Named("R") = R));
}
