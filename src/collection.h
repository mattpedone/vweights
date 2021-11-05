#ifndef __COLLEZIONE_DI_FUNZIONI_HPP__
#define __COLLEZIONE_DI_FUNZIONI_HPP__


#include <RcppArmadillo.h>

#include <cmath>
#include <vector>
#include <string>
/// Per avere infiniro
#include <limits>
#include <iostream>




class COLLEZIONE{
	public:
		//Costructor
		COLLEZIONE(void);
		COLLEZIONE(long double const Loc_kappa, long double const Loc_sig, long double const Loc_omega,  unsigned long int const Loc_n);
		//COLLEZIONE(long double const Loc_kappa0, long double const Loc_sig0, long double const Loc_omega0, long double const Loc_kappa, long double const Loc_sig, long double const Loc_omega,  unsigned long int const Loc_n, unsigned long int const Loc_d);
		
		/* Dovrei fare i setter, ma non li faro'*/

		/*getter Anche questi non li faccio */

		arma::vec main_calcola_Vnk(void);
		arma::vec main_calcola_prior_nclust(void);

		private:
		long double omega;
		long double sig;
		long double kappa;
		unsigned long int n;

};


#endif
