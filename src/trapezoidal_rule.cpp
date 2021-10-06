#include "trapezoidal_rule.h"

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

// #include <gsl/gsl_errno.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_roots.h>
// #include <gsl/gsl_sf_gamma.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_cdf.h>
// 


//// Questa procedura qui sotto non so a cosa serva.
#pragma STDC FENV_ACCESS on



/////////////////////////////// Funzioni per l'integrazione con la regola del Trapezio!

//Costruttore

Integrale_trapezi::
	Integrale_trapezi(long double const &kappa, long double const &sig, long double const &omega, unsigned long int const &n,bool trapezi) 
{
// Parametri NGG:
	M_kappa = kappa;
	M_sig = sig;
	M_omega = omega;
	M_n = n;	


	//std::cout<<__LINE__<<" <- sono a questa linea del file itegrale_con_trapezi.cpp, Ciao"<<std::endl;
	
	//std::cout<<"kappa, sigma, omega, n, t="<<M_kappa<<" "<<M_sig<<" "<<M_omega<<" "<<M_n<<" "<<trapezi<<std::endl;


	int lgr=50000;
	lnintegrale.assign(M_n,0);
	
	long double epsu=0.9999;
	
	//long double app=pow(epsu,(1/(M_n-1)));
	//printf(" vediamo %e\n",epsu-pow(epsu, M_n/(M_n-1)) );
	//std::cout<<"Poi app= "<<log(1-std::pow(epsu,(1/(M_n-1))))<<std::endl;


	long double umax= std::exp( M_sig*std::log(M_omega)-M_sig*(std::log(taylor_f(epsu,15))));
	//std::cout<<"umax="<<log(1-std::pow(epsu,(1/(M_n-1))))<<std::endl;


	if(trapezi){
		for(int kk=1;kk<=M_n;kk++){
			//std::cout<<__LINE__<<" <- sono a questa linea, Ciao"<<std::endl;

			lnintegrale[kk-1]=regola_trapezi(kk,5000,umax);
			//std::cout<<__LINE__<<" <- sono a questa linea, Ciao"<<std::endl;
			//std::cout<<"trapezi["<<kk<<"]="<<lnintegrale[kk-1]<<std::endl;
		}
	}else{//Simpson
		for(int kk=1;kk<=M_n;kk++){
			//std::cout<<__LINE__<<" <- sono a questa linea, Ciao"<<std::endl;

			lnintegrale[kk-1]=simpson(kk,50000,umax);
			//std::cout<<__LINE__<<" <- sono a questa linea, Ciao"<<std::endl;
			//std::cout<<"trapezi["<<kk<<"]="<<lnintegrale[kk-1]<<std::endl;
		}
	}




}

long double 
Integrale_trapezi::
funzione(long double u, int k, bool log, long double logmax){
	
	double out;
	if(M_n==1){
		out=0;
	}
	else{
		out=(M_n-1)*std::log(pow(u,(1/M_sig))-M_omega)-(M_n-1)/M_sig*std::log(u);
	}


	out+=(k-1)*std::log(u)-M_kappa/M_sig*u-logmax;

	//std::cout<<__LINE__<<" <- sono a questa linea, (u,out)="<<u<<","<<out<<std::endl;


	if(log){
		return(out);
	}
	else{
		return(std::exp(out));
	}
}



long double
Integrale_trapezi::
taylor_f(double x,int h){

	double out=(double) (1-x)/(M_n-1);
	//std::cout<<"out= "<< out<<std::endl;
	
	double lout=0;	
	double la=0;
	double lb=0;

	for(int i=2;i<=(h-1);i++){
		la += std::log(i);
		lb += std::log((M_n-1)*i-M_n);
		lout = lb-la+i*std::log(1-x)-i*std::log(M_n-1);
		//std::cout<<"i="<<i<<" lout= "<< lout<<std::endl;
		out += exp(lout);
		//std::cout<<"out= "<< out<<std::endl;

	}

	return(out);
}


long double
Integrale_trapezi::
regola_trapezi(int k, int lgr, double umax, int quante){
	
	
	long double moda = M_sig/M_kappa*(k-1);
	long double lnmax = 0;

	if(moda<pow(M_omega,M_sig)){
		moda=pow(M_omega,M_sig);
	}
	lnmax = (k-1)*log(moda)-M_kappa/M_sig*moda;
	//std::cout<<"lnmax = "<<lnmax<<std::endl;


	double media = k*M_sig/M_kappa;
	double sd = std::sqrt(k)*M_sig/M_kappa;
	//std::cout<<"quante="<<quante<<" media="<<media<<" sd="<<sd<<std::endl;


	double low = media-quante*sd;
	low = std::max(low,pow(M_omega,M_sig));

	double up = media+quante*sd;
	
	//std::cout<<"Prima up= "<<up<<std::endl;


//       f(x):= -ln(1 - x)
//       taylor(f,x);
//                  1  2   1  3   1  4   1  5    / 6\
//              x + - x  + - x  + - x  + - x  + O\x /
//                  2      3      4      5           
//

 	//double umax= std::exp( M_sig*std::log(M_omega)-M_sig*(app+0.5*pow(app,2)+0.333*pow(app,3) ) );
	
	//std::cout<<"Ancora umax= "<<umax<<std::endl;

	up=std::max(up,umax);

	up=std::max(up,3.);
	//std::cout<<"Alla fine up= "<<up<<std::endl;

	/// ORA SONO PRONTO A FARE LA GRIGLIA
	double epsg=(up-low)/(double)lgr;

	//std::cout<<__LINE__<<" <- sono a questa linea, Ciao epsg="<<epsg<<" lgr="<<lgr<<std::endl;

	//std::cout<<"up="<<(up-low)<<" lgr="<<lgr<<std::endl;

	// Applico la formula dei trapezi Vedi Josef Stoer (1974) pg 106 Section 3.2
	double u0=low;
	double integrale=0;
	for(int i=0; i<lgr;i++){
		integrale += funzione(u0,k,false,lnmax);
		u0 += epsg;
	}
	integrale -= (funzione(low,k,false,lnmax)+funzione(up,k,false,lnmax))/2.0;
	integrale *=epsg;
/// Fine formula trapezi

	//std::cout<<"integrale="<<integrale<<std::endl;
	double out = lnmax+std::log(integrale);

	return(out);
}











long double
Integrale_trapezi::
simpson(int k, int lgr, double umax, int quante){


	if ( lgr % 2 != 0 )
    {
	    //std::cout << lgr << " is odd so I will add 1 "<<std::endl;
	    lgr=lgr+1;

    }
	
	long double moda = M_sig/M_kappa*(k-1);
	long double lnmax = 0;

	if(moda<pow(M_omega,M_sig)){
		moda=pow(M_omega,M_sig);
	}
	lnmax = (k-1)*log(moda)-M_kappa/M_sig*moda;
	//std::cout<<"lnmax = "<<lnmax<<std::endl;


	double media = k*M_sig/M_kappa;
	double sd = std::sqrt(k)*M_sig/M_kappa;
	//std::cout<<"quante="<<quante<<" media="<<media<<" sd="<<sd<<std::endl;


	double low = media-quante*sd;
	low = std::max(low,pow(M_omega,M_sig));

	double up = media+quante*sd;
	
	//std::cout<<"Prima up= "<<up<<std::endl;


	
	//std::cout<<"Ancora umax= "<<umax<<std::endl;

	if(moda==1){
	//std::cout<<"moda= "<<moda<<std::endl;
	up=std::max(up,umax);
	}

	up=std::max(up,3.);
	//std::cout<<"Alla fine up= "<<up<<std::endl;

	/// ORA SONO PRONTO A FARE LA GRIGLIA
	double epsg=(up-low)/(double)(lgr-1);

	//std::cout<<__LINE__<<" <- sono a questa linea, Ciao epsg="<<epsg<<" lgr="<<lgr<<std::endl;

	//std::cout<<"up="<<(up-low)<<" lgr="<<lgr<<std::endl;

	// Applico la formula di  Simpson. L'ho presa su Wikipedia... :-/ 
	double u1,u2,u3;
	double integrale=0;


	for(int j=1; j<=(lgr/2);j++){
		u1=low+(2*j-2)*epsg;
		u2=low+(2*j-1)*epsg;
		u3=low+(2*j-0)*epsg;
		integrale += epsg/3*(funzione(u1,k,false,lnmax)+4*funzione(u2,k,false,lnmax)+funzione(u3,k,false,lnmax) ) ;

		if(std::isnan(integrale)){ 
			//std::cout<<"j="<<j<<" u1="<<u1<<" u2="<<u2<<" u3="<<u3 << " integrale="<<integrale<<std::endl;
			/*std::cout<<std::log(pow(u1,(1/M_sig))-M_omega)+0.<<std::endl;
			std::cout<< (M_n-1)/M_sig*std::log(u1)<<std::endl;
			std::cout<<(M_n-1)*std::log(pow(u1,(1/M_sig))-M_omega)-(M_n-1)/M_sig*std::log(u1)<<"\n";
			std::cout<<(k-1)*std::log(u1)-M_kappa/M_sig*u1-lnmax<<"\n";*/
		
		
		}
		}
	
	/// Fine formula Simpson

	//std::cout<<" integrale="<<integrale<<std::endl;

	double out = lnmax+std::log(integrale);

	return(out);
}






