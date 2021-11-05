#include "collection.h"
#include "trapezoidal_rule.h"
#include "genstirling_rec.h"

#include <RcppArmadillo.h>

#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <random>
#include <iostream>
#include <cfenv>
#include <algorithm>
/// Per avere infiniro
#include <limits>


#include <gsl/gsl_sf_gamma.h>

//// Questa procedura qui sotto non so a cosa serva.
#pragma STDC FENV_ACCESS on


COLLEZIONE:: COLLEZIONE(void){
	//std::cout<<"Hai chiamato il costruttore di COLLEZIONE vuoto"<<std::endl;
}

COLLEZIONE::
COLLEZIONE(long double const Loc_kappa, long double const Loc_sig, long double const Loc_omega,  unsigned long int const Loc_n){

	//Parametri NGG:
	kappa = Loc_kappa;
	sig= Loc_sig;
	omega= Loc_omega;

	//data
  n = Loc_n;
}


arma::vec
	COLLEZIONE::main_calcola_Vnk(void){
		{

		//This line should be all, what you need to use infinity in the code.
	  static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
			// Here I compute the log(integral) in the V_{k,n} formula by a
			// Simpson quadrature rule!
			Integrale_trapezi logsimpson(kappa,sig,omega,n);
			
			arma::vec Vnk(n);
			
			for(int kk=1;kk<=n;kk++){
			  Vnk[kk-1]=std::exp(kappa/sig*std::pow(omega,sig)+kk*std::log(kappa)-(+1)*std::log(sig)-gsl_sf_lngamma(n/1.0)
                          +logsimpson.output()[kk-1]);
			  //std::cout<<Vnk[kk-1] <<std::endl;
		}
			return Vnk;
	}
	  
	}

arma::vec
  COLLEZIONE::main_calcola_prior_nclust(void){
    {
    //This line should be all, what you need to use infinity in the code.
    static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
   //This class uses a recursive formula to compute the generalized
    //Stirling number of second kind in log scale.
    Gen_stirling logStirling(n,sig,true);

    //The vector of Stirling number S(n,k,sigma) k=1,...,n
    //is saved in a vector via a getter, be careful because
    // vectors in C++ start from 0.
    
    std::vector<long double> logSt=logStirling.output();
    
    // Here I compute the log(integral) in the V_{k,n} formula by a
    // Simpson quadrature rule!
    Integrale_trapezi logsimpson(kappa,sig,omega,n);
    // Questo dovrebbe aprire un file di uscita nella catella path
    //std::ofstream FILE(path, std::ios::out | std::ofstream::binary);
    
    //long double cumulata=0;
    //std::vector<long double> prob(n,0);
    arma::vec prob(n);
    for(int kk=1;kk<=n;kk++){
      prob[kk-1]=std::exp(kappa/sig*std::pow(omega,sig)+kk*std::log(kappa)-(kk+1)*std::log(sig)-gsl_sf_lngamma(n)+logsimpson.output()[kk-1]+logSt[kk-1]);
      //cumulata += prob[kk-1];
      //std::cout<<prob[kk-1] <<std::endl;
      }
    
    //std::cout<<"La cumulata è pari a "<< cumulata <<std::endl;

    // I write the vector of long double prob in the file prob.out
    //UTIL::Ldouble2file(prob,"prob.out");
   
    return prob;
    }
  }
