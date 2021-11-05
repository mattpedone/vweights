#include <RcppArmadillo.h>

#include <random>
#include <functional>
#include <iostream>
#include <fstream>
#include <cmath>

#include <sstream>
#include <gsl/gsl_sf_gamma.h>

#include "trapezoidal_rule.h"
#include "gensterling_rec.h"
#include "collection.h"

#include <vector>
#include <string>

//' Computes the weights
//'
//' @param Locn Number of observation
//' @param Locsig Sigma parameter
//' @param Lockappa Kappa parameter
//' @return vector of length Locn
//' @export
// [[Rcpp::export]]
arma::vec computev(int Locn, double Locsig, double Lockappa){
  //This function calls a main! Here I only set the parameters

  ///////////////

  //Locn    // The sample size
  //Locsig  // The sigma parameter
  //Lockappa// the kappa parameter

  ///////////////

  double Locomega = 1.0; // Actually omega is a scale parameter. See the
  // Argiento Guglielmi Pievatolo (2010) CSDA, we can fix it at 1

  COLLEZIONE collezione(Lockappa,Locsig,Locomega,Locn);

  arma::vec Vnk(Locn);
  Vnk = collezione.main_calcola_Vnk();

  return Vnk;
}

//' Computes the prior probability mass function for the number of clusters
//'
//' @param Locn Number of observation
//' @param Locsig Sigma parameter
//' @param Lockappa Kappa parameter
//' @return vector of length Locn
//' @export
// [[Rcpp::export]]
arma::vec compute_pnclu(int Locn, double Locsig, double Lockappa){
  //This function calls a main! Here I only set the parameters
  
  ///////////////
  
  //Locn    // The sample size
  //Locsig  // The sigma parameter
  //Lockappa// the kappa parameter
  
  ///////////////
  
  double Locomega = 1.0; // Actually omega is a scale parameter. See the
  // Argiento Guglielmi Pievatolo (2010) CSDA, we can fix it at 1
  
  COLLEZIONE collezione(Lockappa,Locsig,Locomega,Locn);
  
  arma::vec prob(Locn);
  prob = collezione.main_calcola_prior_nclust();
  
  return prob;
}