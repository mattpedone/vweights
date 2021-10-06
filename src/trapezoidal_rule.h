#ifndef __INTEGRALE_CON_TRAPEZI_HPP__
#define __INTEGRALE_CON_TRAPEZI_HPP__



#include <cmath>
#include <vector>
#include <string>
/// Per avere infiniro
#include <limits>
 

#include <iostream>

class Integrale_trapezi{

	public:
		//Costruttore
		Integrale_trapezi(){
			//std::cout<<"Ciao dal costruttore Trapezi vuoto"<<std::endl;
		};
		Integrale_trapezi(long double const &kappa, long double const &sig, long double const &omega, unsigned long int const &n,bool trapezi=false); 
 		~Integrale_trapezi(){
 			//std::cout << "Distruggo la classe trapezi! Non capisco se ha senso..." << std::endl;
 	};

		/* getters */
		std::vector<long double> output(){ return lnintegrale;}

		//*  Funzione di utilità * mi sa che la posso mettere pubblica
		long double regola_trapezi(int k, int lgr,double umax, int quante=40);
		long double simpson(int k, int lgr, double umax, int quante=40);

		// //Compute the taylor expansion around 1 of (1-x^(n-1))
                long double taylor_f(double x,int h);


	private: 
		// Parametri del NGG:
		long double M_omega;
		long double M_sig;
		long double M_kappa;
		unsigned long int M_n;   // numerosità campione
		
		
		/// Variabile di uscita
		std::vector<long double> lnintegrale;
		
		
		// Funzioni!
		long double funzione(long double u, int k, bool log, long double logmax=0);
};//Chiudo la classe




#endif
