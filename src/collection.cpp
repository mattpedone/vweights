#include "collection.h"
#include "trapezoidal_rule.h"

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
