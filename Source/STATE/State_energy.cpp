#include "State.h"

using namespace std;
void State::CalE0()
{
  E0 = 0;
  
  if(!network) {
    for(int i=0; i<N; i++)  {
      GetNeighbors(i);
      
      int nsame=0;
      for(int j=0; j<Z; j++)
	if(spin[i] == spin[neighborLocs[j]]) 
	  nsame++;
      
      E0 += D-nsame - h[i]*spin[i];

    }
  } 
  

  
  else { 
    for(int i=0; i<N; i++)  {

      E0 -= h[i]* spin[i];
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = *p;
	double Jij = edgeweight(i, j); 

	E0 -= 0.5* Jij * spin[i] * spin[j];
      }
    }

  }
}

