/*######################################################################################################################
  #  Generate synthetic node weights and edge weights to test the RFIM and MWCS algorithms.                              #
  #  Yang-Yu Liu                                                                                                            #
  #  02/10/15                                                                                                            #
  #######################################################################################################################*/
#include <iostream>
#include <list>
#include <queue>
#include "types.h"   
#include "hi_pr.h" 
#include "State.h"   
#include "LineSegment.h"
#include "CrossingField.h"               
#include "DropletAnalyse.h"    
#include "MakeDirectory.h"

using namespace std;
//////////////////////////////////      GLOBAL CONSTANT          ///////////////////////////////////////////////////////
//const hType  RESOLUTION=(hType)1e15;// used to distinguish two arbitrary random fields in the lattice
const hType  RESOLUTION=(hType)1e17;// used to distinguish two arbitrary random fields in the lattice
                                    // resolution = 1/discreteness

int dist;                       // the random field distribution
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int main (int argc, char** argv)
{
  if(argc!=5)
    {
      cout <<"\nNumber of arguments is incorrect.";
      cout <<"\nGenerate synethic node weights and edge weights on a random network constructed with static model.";
      cout <<"\nCommand format:  ./HMS_GenerateSyntheticData N c gamma seed\n";
      exit(0);
    }

  float t0   = timer();                   // get the beginning time    

  int N        = atoi(argv[1]);             // # of nodes
  double c     = atof(argv[2]);             // mean degree 
  double gamma = atof(argv[3]);             // degree exponent used in the static model
  int seed     = atoi(argv[4]);             // random seed

  
  cout<<"-------------------------------------------------------------------------------------------------------\n";
  cout <<"Generate synethic node weights and edge weights on a random network constructed with static model.\n";
  cout << " N="<< N <<" c="<< c << " gamma= " << gamma << " seed=" << seed << endl;
  MakeDirectory_Default();


  State C(N);  
  C.Set_RandomPvalue_RandomBond(N, c, gamma, seed);
  
  float t = timer() - t0;                
  cout << "\n It takes " << t << " s \n";

}
