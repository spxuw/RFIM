#include "State.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////
/*  Attention:
    When we calculate the energy of each spin site, we should consider the energy of its Z 
    bonds actually is shared by its Z neighbors. 
    Thus, 
 
    Ei = 0.5*[-J * (nsame*1 + (Z-nsame)*(-1))]
       = 0.5*[J*(Z-2 nsame)]
       = J*(D-nsame)

    For network with random bonds,
    Ei = 0.5*[-Jij * Si*Sj]

    Note that the total energy 
    E = - sum_{i~j} Jij Si Sj - sum_{i} hi Si 

     
*/
void State::CalE0()// just used to calculate the internal energy of the GS
{
  E0 = 0;   //E=0;    //M=0;
  
  if(!network) {
    //cout << "For lattice " << endl; //debug
    for(int i=0; i<N; i++)  {
      GetNeighbors(i);
      
      int nsame=0;  //count how many neighbors have the same spin
      for(int j=0; j<Z; j++)
	if(spin[i] == spin[neighborLocs[j]]) 
	  nsame++;
      
      E0 += D-nsame - h[i]*spin[i];          // if spin[i]= -1 or +1
      //E0 += D-nsame - h[i]*(2*spin[i]-1);  // if spin[i]=  0 or +1

    }
  } // end of if lattice
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  else { // for network
    //cout << "For network " << endl; //debug

    for(int i=0; i<N; i++)  {

      E0 -= h[i]* spin[i];
      //E0 -= h[i]* (2*spin[i]-1); 
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = *p;
	double Jij = edgeweight(i, j); 
	//cout << Jij << '#'; //test

	E0 -= 0.5* Jij * spin[i] * spin[j];                  // if spin[i]= -1 or +1
	//E0 -= 0.5* Jij *  (2*spin[i]-1) *  (2*spin[j]-1);  // if spin[i]=  0 or +1
      }
    }

  }// end of if network 
  //cout << "E0 = " << E0 << endl; //test 

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
