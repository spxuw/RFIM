#include <iostream>
#include <fstream>
#include <vector>

#include "Rand.h"
#include "gsea.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if(argc!=4) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : test_GSEA N M seed\n ";
    exit(0);
  }
   
  int N    = atoi(argv[1]);
  int M    = atoi(argv[2]);
  int seed = atoi(argv[3]);

  Rand rand;
  rand.seed(seed);

  vector<double> r0(N);
  for(int i=0; i<N; i++) 
    r0[i] = 2*rand.ran1()-1; // randomly generate a correlation coefficient [-1,1]

  vector<int> S0;
  //rand.randomsubset_knuth(N, M, S0);

  double p = M/(double)N;
  for(int i=0; i<N; i++) {
    if(i<0.2*N) {
      if(rand.ran1()<5*p)
	S0.push_back(i);
    }
    else if (i<0.7*N) {
      if(rand.ran1()<p)
	S0.push_back(i);
    }
    else {
      if(rand.ran1()<0.2*p)
	S0.push_back(i);
    }
  }
  
  M = S0.size();
  for(int j=0; j<M; j++) 
    cout << S0[j] << endl;


  CalES(r0, S0);

  exit(1);
}

