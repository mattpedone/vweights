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
    
    int lgr=50000;
    lnintegrale.assign(M_n,0);
    
    long double epsu=0.9999;
    
    long double umax= std::exp( M_sig*std::log(M_omega)-M_sig*(std::log(taylor_f(epsu,15))));
    
    if(trapezi){
      for(int kk=1;kk<=M_n;kk++){
        
        lnintegrale[kk-1]=regola_trapezi(kk,5000,umax);
      }
    }else{//Simpson
      for(int kk=1;kk<=M_n;kk++){
        lnintegrale[kk-1]=simpson(kk,50000,umax);
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
      
      double lout=0;	
      double la=0;
      double lb=0;
      
      for(int i=2;i<=(h-1);i++){
        la += std::log(i);
        lb += std::log((M_n-1)*i-M_n);
        lout = lb-la+i*std::log(1-x)-i*std::log(M_n-1);
        out += exp(lout);}
      
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
      
      double media = k*M_sig/M_kappa;
      double sd = std::sqrt(k)*M_sig/M_kappa;
      
      double low = media-quante*sd;
      low = std::max(low,pow(M_omega,M_sig));
      
      double up = media+quante*sd;
      
      //       f(x):= -ln(1 - x)
      //       taylor(f,x);
      //                  1  2   1  3   1  4   1  5    / 6\
      //              x + - x  + - x  + - x  + - x  + O\x /
      //                  2      3      4      5           
      //
      
      up=std::max(up,umax);
      
      up=std::max(up,3.);
      
      /// ORA SONO PRONTO A FARE LA GRIGLIA
      double epsg=(up-low)/(double)lgr;
      
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
      
      double out = lnmax+std::log(integrale);
      
      return(out);
    }

long double
  Integrale_trapezi::
    simpson(int k, int lgr, double umax, int quante){
      
      
      if ( lgr % 2 != 0 )
      {
        lgr=lgr+1;
      }
      
      long double moda = M_sig/M_kappa*(k-1);
      long double lnmax = 0;
      
      if(moda<pow(M_omega,M_sig)){
        moda=pow(M_omega,M_sig);
      }
      lnmax = (k-1)*log(moda)-M_kappa/M_sig*moda;
      
      double media = k*M_sig/M_kappa;
      double sd = std::sqrt(k)*M_sig/M_kappa;
      
      double low = media-quante*sd;
      low = std::max(low,pow(M_omega,M_sig));
      
      double up = media+quante*sd;
      
      if(moda==1){
        up=std::max(up,umax);
      }
      
      up=std::max(up,3.);
      
      /// ORA SONO PRONTO A FARE LA GRIGLIA
      double epsg=(up-low)/(double)(lgr-1);
      
      // Applico la formula di  Simpson. L'ho presa su Wikipedia... :-/ 
      double u1,u2,u3;
      double integrale=0;
      
      for(int j=1; j<=(lgr/2);j++){
        u1=low+(2*j-2)*epsg;
        u2=low+(2*j-1)*epsg;
        u3=low+(2*j-0)*epsg;
        integrale += epsg/3*(funzione(u1,k,false,lnmax)+4*funzione(u2,k,false,lnmax)+funzione(u3,k,false,lnmax));
      }
      
      /// Fine formula Simpson
      
      double out = lnmax+std::log(integrale);
      
      return(out);
    }
