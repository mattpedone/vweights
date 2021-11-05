#ifndef __GEN_STIRLING_HPP__
#define __GEN_STIRLING_HPP__

#include <cmath>
#include <vector>
#include <string>
/// Per avere infiniro
#include <limits>

#include <iostream>

// Funzioni per il calcolo dell generalized stirling.
// Voglio utilizzare la funzione ricorsiva di Charalambides
// cp.2 pg 104 Th 218 formula 2.67
//
// Dove per noi r=0 e se S(n,k,sig) è il generalized stirling di "Controlling..." allor
// dobbiamo pensare che S(n,k,sig)=(-1)^(n-k)C(n,k,sig,0)


// In particolare mi creerei il vettore di lunghezza n che conitiene S(n,k,sig) k=0,...,n
//
class Gen_stirling
{
public:
  //Costruttore
  Gen_stirling(){
    std::cout<<"Costruttore vuoto Stirling"<<std::endl;
  }
  Gen_stirling(unsigned long int const &n, long  double const &sig);
  // Costruttore specificando log scale
  Gen_stirling(unsigned long int const &n, long  double const &sig,bool const &log_scale);
  ~Gen_stirling(){
    ////std::cout<<"Ciao dal DISTRUTTORE di Stirling"<<std::endl;
  }
  
  //this function apply  formula 267 modificata in log scale to the vector da aggiornare 
  //
  void fai_un_passo(std::vector<long double> &da_aggiornare);
  /* getters */
  std::vector<long double> output(){ 
    std::vector<long double>::const_iterator first = precedente.begin() + 1;
    std::vector<long double>::const_iterator last = precedente.end();
    std::vector<long double> uscita(first,last);
    //std::cout<<"Last è "<<*last <<std::endl;
    return uscita;
    
  }
  
  /* printer non so se si usa */
  void a_che_passo_sono(){
    std::cout<<"Sono al passo="<<passo<<std::endl;
  }
  
private: 
  unsigned long int  M_n;  // Primo argomento dello Gen Stirling
  long double  M_sig;      // Terzo argomento dello Gen Stirling
  
  bool M_flag;// This variable is needed to understand if we want the log or the natural scale of the Stirling
  double out;
  
  std::vector<long double> precedente;
  std::vector<long double> successivo;
  
  int passo;
  
  // This function implement formula 267 of Charambides modified 
  // by the (-1)^(n-k) factor so it holds for the S(n,k,sig) representation
  void formula267_modificata(void);  
  void calcola_ultimo();		
  
  // This function implement formula 267 in log scale Charambides modified 
  // by the (-1)^(n-k) factor so it holds for the S(n,k,sig) representation
  void formula267_log(void); 
};

// In particolare mi creerei il vettore di lunghezza n che conitiene S(n,k,sig) k=0,...,n

//Definisco una FUNZIONE! che mi torna, per referenza, 
//la matrice dei numeri di Stirling (M_n,M_n)
// I know that i C++ function are  not the best thingh to do. 
// But in this case is very naural

void
  Gen_lstirling_matrice(unsigned long int const &n, long  double const &sig, std::vector<std::vector<long double>> &LogSirling);



/////////////// gsl function: For the integration  //////////////////
#include <gsl/gsl_integration.h> 
struct f_params { long double kappa; long double sig; long double omega; unsigned long int n; unsigned long int k; double lnmax;};

double integranda (double u, void *p);  // p sono i parametri, u la variabile della funzione
/////////////////////////////////////////////////////////////////////////////////////////////////////

class Gen_integrale
{
public:
  //Costruttore
  Gen_integrale(long double const &kappa, long double const &sig, long double const &omega, unsigned long int const &n); 
  Gen_integrale(long double const &kappa, long double const &sig, long double const &omega, unsigned long int const &n,bool const &ciclica);	
  /* getters */
  inline std::vector<double> output()  { return lnintegrale; }
  
  //*  Funzione di utilità * mi sa che la posso mettere pubblica
  long double calcola_lnintegrale(unsigned long int k,unsigned long int nn);
  
private: 
  // Parametri del NGG:
  long double M_omega;
  long double M_sig;
  long double M_kappa;
  unsigned long int M_n;   // numerosità campione
  long double Gamma;       // Salvo la Gamma di n
  
  /// Variabile di uscita
  std::vector<double> lnintegrale;
  // Variabile che mi serve solo quando uso la forma ciclica
  std::vector<double> lnsuccessivo;
  /// varialbile che mi serve per fare il ciclo 
  int passo;
  
  // Variabili per la formula di integrazione:
  double epsabs;  // errore assoluto massimo
  double epsrel;  // errore relativo massimo
  size_t limit;   // numero massimo di intervalli usati dalla procedura adattiva
  size_t NI;      // numero di intervalli riservati in memoria nella formula di quadratura
  int key; 	// Vedere http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html#QAG-adaptive-integration
  // per i nodi scelti 
  gsl_integration_workspace * WS; // Crea una griglia di nodi
  
  // Funzioni per la forma ciclica
  void formula_mia(void);
  void cicla_Vnk(void);
};

#endif
