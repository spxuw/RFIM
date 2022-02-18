#include "State.h"

using namespace std;



void State::HMS_Ideker(int itype, double Ti, double Tf, int seed, char* file)  
{
  cout << "\nSearch for the high-scoring modules using Ideker's simulated annealing algorithm: " << endl;

  char fname [256];
  sprintf(fname,"%s.Score_T", file);
  ofstream fout(fname, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  rand.seed(seed); 
  
  if(itype==1) {
    for(int i=0; i<N; i++) 
      spin[i] = (rand.ran1()<0.5) ? ACTIVE : INACTIVE;
  }
  else if(itype==2) {
    for(int i=0; i<N; i++) 
      spin[i] = (rand.ran1()>PVALUE[i]) ? ACTIVE : INACTIVE; 
  }
 


  
  
  

  double S = GetScore();

  int iteration = 0;
  for(double T=Ti; T>=Tf; T*=0.5, iteration++) { 
    cout << "T=" << T << ", S=" << S << endl;
    fout << iteration << ' ' << T << ' ' << S << endl;

    for(int t=0; t<N; t++) {
      int i = rand.discrete(0,N-1);
      spin[i] *= (-1);

      
      
      

      double Snew = GetScore();
  
      if(Snew>S) { 
	S = Snew;
      }
      else { 
	if(rand.ran1() < exp((Snew-S)/T))
	  S = Snew; 
	else 
	  spin[i] *= (-1);
      }
    } 

  }
  fout.close();
 
  
  
}


