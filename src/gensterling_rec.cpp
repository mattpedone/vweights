#include "genstirling_rec.h"

#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include <random>
#include <iostream>
#include <cfenv>
#include <algorithm>
/// Per avere infinito
#include <limits>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


//// Questa procedura qui sotto non so a cosa serva.
#pragma STDC FENV_ACCESS on

// Definisco ora il costruttore a 3 parametri, in questo caso
// posso scegliere la scala logaritmica
Gen_stirling::
  Gen_stirling(unsigned long int const &n, long  double const &sig, bool const &log_scale){
    M_n= n;
    M_sig=sig;
    
    M_flag=log_scale; // In questo caso faccio il conto nella scala naturale
    
    if(!M_flag){
      Gen_stirling App(M_n,M_sig);
      precedente=App.output();
      return ;
    }
    
    double infinito = std::numeric_limits<double>::infinity();
    // Initialize passo
    passo=1;
      
      // Initialize precedete and successivo
      precedente.assign(M_n+1,-infinito); // Quando faccio il conto in scala naturale
      precedente[1]=std::log(M_sig);
      //parto dal vettore con tutti zeri tranne il primo elemento che è uno
      
      // Successivo lo inizializzo a meno infinito!
        successivo.assign(M_n+1,-infinito);
        
      //chiama la funzione principale out;
      calcola_ultimo();
}


/// Charalambides 
// This function implement formula 267 of Charambides modified 
// by the (-1)^(n-k) factor so it holds for the S(n,k,sig) representation

void
Gen_stirling::
formula267_modificata(void){

	// In entrata nn è pari a passo
	int nn= passo;
	// ora posso incrementare passo
	passo=passo+1;

	successivo[0]=0;
	for(int kk=1;kk<=passo;kk++){
		successivo[kk]=(-M_sig*kk+nn)*precedente[kk]+M_sig*precedente[kk-1];

	}
      precedente.swap(successivo);
}

// This function implement formula 267 in log scale Charambides modified 
// by the (-1)^(n-k) factor so it holds for the S(n,k,sig) representation
void 
Gen_stirling::
formula267_log(void){
	// In entrata nn è pari a passo
	int nn=passo;
	// ora posso incrementare passo
	passo=passo+1;

	double infinito = std::numeric_limits<double>::infinity();

	successivo[0]=-infinito;
	for(int kk=1;kk<passo;kk++){
		successivo[kk]=std::log(nn-M_sig*kk)+precedente[kk]+
			std::log(1+M_sig/(nn-M_sig*kk)*std::exp(precedente[kk-1]-precedente[kk]));

	}
	successivo[passo]= log(M_sig)+precedente[nn];

	precedente.swap(successivo);
}

void 
Gen_stirling::
fai_un_passo(std::vector<long double> &da_aggiornare){

	// In entrata nn è pari alla dimensione attuale di da_aggiornare
	int nn=da_aggiornare.size();
	// ora posso incrementare nn
	int da_a=nn+1;

	double infinito = std::numeric_limits<double>::infinity();
	std::vector<long double> appoggio(da_aggiornare);

	da_aggiornare.push_back(0.);	
	da_aggiornare[0]=std::log(nn-M_sig*1)+appoggio[0];
	{
		int kk;
		for(int ki=1;ki<(da_a-1);ki++){
			kk=ki+1;
				da_aggiornare[ki]=std::log(nn-M_sig*kk)+appoggio[ki]+
				std::log(1+M_sig/(nn-M_sig*kk)*std::exp(appoggio[ki-1]-appoggio[ki]));

		}
	}
	da_aggiornare[da_a-1]= log(M_sig)+appoggio[nn-1];
}

void
Gen_stirling::
calcola_ultimo(void){

	if(!M_flag){
		while(passo<M_n){
			formula267_modificata();
		}
	}
	if(M_flag){
		while(passo<M_n){
			formula267_log();
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
/////////////////// Now I will define a function 
//Definisco una FUNZIONE! che mi torna, per referenza, 
//la matrice dei numeri di Stirling (M_n,M_n)
// I know that i C++ function are  not the best thingh to do. 
// But in this case is very naural

void
Gen_lstirling_matrice(unsigned long int const &n, long  double const &sig, std::vector<std::vector<long double>> &LogStirling){

	int M_n = n;
	long double M_sig = sig;
	
	std::cout<<__LINE__<<" CIAO DALL FUNZIONE per la Matrice"<<std::endl;
	double infinito = std::numeric_limits<long double>::infinity();

	// Initialize the first row of the matrix LogStirling

	std::cout<<__LINE__<<" CIAO DALL FUNZIONE per la Matrice"<<std::endl;

	std::cout<<(LogStirling[0]).size()<<std::endl;

	std::cout<<__LINE__<<"(LogStirling[0])[0]";
	std::cout <<"\n";
	std::cout << (LogStirling[0])[0]<<std::endl;

	(LogStirling[0])[0] =std::log(M_sig);
	
	std::cout<<__LINE__<<" CIAO DALL FUNZIONE per la Matrice"<<std::endl;
	std::cout << (LogStirling[0])[0]<<std::endl;

	//parto dal vettore con tutti zeri tranne il primo elemento che è uno

	//chiama la funzione principale out;

	int nn,np1,kk;
	for(int riga=1;riga<=(M_n-1);riga++){
		// Inizializzo la successiva riga della matrice !
		np1=riga+1; //questo è n+1 degli appunti 
		nn= riga;    //questo è n degli appunti     

		//LogStirling[riga].assign(-infinito,np1);

		LogStirling[riga][0]= std::log(riga-M_sig)+LogStirling[riga-1][0];
		for(int kc=1;kc<=(riga-1);kc++){
			kk=kc+1;
			// kk è il k delgi appunti

			LogStirling[riga][kc]=std::log(riga-M_sig*kk)+ LogStirling[riga-1][kc]+
				std::log(1+M_sig/(nn-M_sig*kk)*std::exp(LogStirling[riga-1][kc-1]-LogStirling[riga-1][kc]));
		}
		LogStirling[riga][riga]= log(M_sig)+LogStirling[riga-1][riga-1];
	}



for(int i=0;i<M_n;i++){
	for(int j=0;j<M_n;j++){
		std::cout<< LogStirling[i][j]<<" ";
	}
	std::cout<<"\n";
}
}




////////////////////////////// INIZIAMO A DEFINIRE LE FUNZIONI DELLA NUOVA CLASSE 


//Prima pero' definisco la funzione da integrare
// Funzione integranda: serve per costruire una gsl_function
// L'integranda è la funzione
//    (1-omega/(u^(1/sig)))^(n-1)*u^(k-1)*exp(-kappa/sig*u)
// Vedi articolo "Controlling the reinforcement.." di Lijoi et al, Appendix.
// ATT: l'integranda potrebbe essere difficile da valutare quando sigma assume valori molto piccoli!
  double
integranda (double u, void *p) 
{
  f_params* params = static_cast<f_params*>(p);
  unsigned long int n=(params->n);
  unsigned long int k=(params->k);
  long double sig=(params->sig);
  long double kappa=(params->kappa);
  long double omega=(params->omega);
  double lnmax=(params->lnmax);
  double out;
  double pezzettino;
  
  pezzettino = 1./sig*log(u);
  // POSSIBILE APPROSSIMAZIONE: quando u^(1/sig) è maggiore di e^100, 
  // si approssima log(u^(1/sig)-omega) con log(u^(1/sig))...
  // In tal caso l'errore introdotto sarebbe:
// -sum_k=(1...+inf){1/(k*u^(k/sig))}   
// (basta vedere sezione Logritmo del librone degli integrali, ln(x/(1-x))) 

//	if(pezzettino > 100) 
//	{
//		out = (k-1)*log(u) - kappa/sig*u - lnmax;
//	}
//	else
	double out1;
	{
		//std::cout<<"frazione"<<1-exp(-pezzettino)<<std::endl;
		out1 = (n-1)*log(1-omega/(pow(u,1./sig)))+(k-1)*log(u)-kappa/sig*u;
	}
	
	//std::cout<<"out1="<<out1<<std::endl;
	out1=out1-lnmax;
	out=exp(out1);
	//std::cout<<"out="<<out<<std::endl;

	if(out>1E301){
		out=1E300;
		std::cout<<"ATT: +inf and u="<<u<<" k="<<k<<" lnmax="<<lnmax<<" out1="<<(k-1)*log(u)-kappa/sig*u<<std::endl;
	}
	if(out<-1E301){
		out=-1E300;
		std::cout<<"ATT: -inf"<<std::endl;
	
	}

	return out;
	
};

// Devo definire il costruttore della nuova clasee
Gen_integrale::
	Gen_integrale(long double const &kappa, long double const &sig, long double const &omega, unsigned long int const &n) 
{
// Parametri NGG:
	M_kappa = kappa;
	M_sig = sig;
	M_omega = omega;
	M_n = n;	
 
// Variabili per la formula di integrazione:
	epsabs= pow(10, -5);  // errore assoluto massimo
	epsrel= pow(10, -5);  // errore relativo massimo
	limit=1e5;        // numero massimo di intervalli usati dalla procedura adattiva
	NI=2e5;           // numero di intervalli riservati in memoria nella formula di quadratura
	key=4; 	              // Vedere http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html#QAG-adaptive-       integration per i nodi scelti: può stare tra 1 e 6
	WS=gsl_integration_workspace_alloc (NI); // Crea una griglia di nodi
	
	lnintegrale.assign(M_n,0);

	for(int kk=1;kk<=M_n;kk++){
		lnintegrale[kk-1]=calcola_lnintegrale(kk,M_n);
		//std::cout<<"Vnk["<<M_n<<","<<kk<<"]="<<M_kappa/M_sig*std::pow(M_omega,M_sig)+kk*std::log(M_kappa)-std::log(M_sig)-gsl_sf_lngamma(M_n+1)+lnintegrale[kk-1]<<std::endl;
		// std::cout<<"lnintegrale["<<kk<<"]="<<lnintegrale[kk-1]<<std::endl;
		//std::cout<<"In formula normalizzazion è "<<M_kappa/M_sig*std::pow(M_omega,M_sig)+kk*std::log(M_kappa)-std::log(M_sig)-gsl_sf_lngamma(M_n+1)<<std::endl;
	}


	// Penso che va questa cosa ma forse posso proprio eliminare WS
	gsl_integration_workspace_free(WS); 
}

long double
Gen_integrale::
calcola_lnintegrale(long unsigned int k, long unsigned int nn){
	
	
	gsl_integration_workspace * w_gsl= gsl_integration_workspace_alloc(9000);
	// Integrazione:
	double integrale(0.0);
	double ERR(0.0);

	gsl_function  F;
	double moda = M_sig/M_kappa*(k-1);
	double lnmax = 0;

	double max=(k-1)*M_sig/M_kappa;
	if(max<pow(M_omega,M_sig)){
	lnmax = -M_kappa/M_sig;
	}
	else{
	lnmax = (k-1)*log(max)-M_kappa/M_sig*max;
	//std::cout<<"lnmax = "<<lnmax<<std::endl;
	}

	// A COSA SERVE lnmax: per stabilizzare il calcolo dell'integrale \int_1...inf f(x)dx, infatti
  //            I = \int_1...inf f(x)dx = M * \int_1...inf {f(x)/M} dx = M * I_s
  //           logI = logM + log(I_s)
  // Scegliamo M = - kappa/sig che è il valore del logaritmo della densità di una gamma(k, kappa/sig) quando u=1,
  // così da fare uno zoom intorno alla distribuzione quando u = 1, cioè sulla coda.
  // In tal modo il comportamento dell'integrale migliora, dato che quando sig è piccolo l'integranda
  // è all'incirca concentrata intorno a sig*k/kappa con una varianza piccolissima (basta considerare che in tal caso 
	// 1-omega/u^(1/sig)) è circa 1, quindi l'integranda è la densità di una gamma appunto. 
  
  /// Proviamo a fare un ragionamento analogo ma pensando gia' all'integrale. Il massimo che puo assumere l'integrale
	//  e dato da sigma^k/kappa^k*Gamma_inc(k,a/sigma*omega^sigma) allora diamo il log di questo valore a lmax
	// Ma dato che la Gamma_incompleta è difficile da calcolare facciamo
	//  Gamma_inc(k,kappa/sigma*omega^sigma)<Gamma(k)=(k-1)!
	
// 	if(k>1000){
// 		std::cout<<"Prima lnmax = "<<lnmax<<std::endl;
// 		lnmax=k*log(M_sig)-k*log(M_kappa)+gsl_sf_lnfact(k-1);
// 		std::cout<<"Dopo lnmax = "<<lnmax<<std::endl;
// 
// 	}
// 

	f_params params = { M_kappa, M_sig, M_omega, nn, k, lnmax }; 

	F.function = &integranda;
	F.params = static_cast<void*>(&params);	// Setto i parametri della funzione

	// Usiamo la funzione della gsl per gli integrali da a a +inf. L'algoritmo è adattivo.
  // gsl_integration_qagiu (gsl function * f, double a, double epsabs, double epsrel, size t limit, gsl integration workspace *
                              // workspace, double * result, double * abserr)
  gsl_integration_qagiu(&F, pow(M_omega,M_sig),  0., 1e-5, 9000, w_gsl, &integrale, &ERR);
  
  double lnintegrale;
  /// Attenzione dato che ho sottratto lnmax al log della funzione integranda
  //è come se avessi moltiplicato per exp(-lnmax) l'integrale
	//quindi devo sommare lnmax in scala logaritmica:
	//std::cout<<"Log Integrale = "<<std::log(integrale)<<std::endl;

	lnintegrale = lnmax + std::log(integrale);

		if(ERR>std::pow(10,-2))
	{
	std::cout<<"Warning: the error of the integration is greater than "<<std::pow(10,-2)<<": ERR= "<<ERR<<
		" --------> lnintegrale("<<k<<")=  "<<lnintegrale<<std::endl;
	}
	gsl_integration_workspace_free(w_gsl);
	return lnintegrale;
}

