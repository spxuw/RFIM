#include "State.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////// 
// default constructor
State::State() 
{
  network=false;

    D = 1;
    L = 1;
    Z = 2;
    N = 1;
    R = 1.0;
    seed = 1;
    Hext = 0;
    //E0 = 0;
    //E  = 0;
    M  = N;
    m  = 1.0;

    spin = new char [N];
    nUp = new int [N];
    for(int i=0;i<N;i++) {
	spin[i] = DOWN;
	nUp[i] = 0;
    } 
    
    avalancheindex = new int [N];

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    for(int i=0;i<Z;i++)
	neighborCoords[i] = new int[D];

    for(int i=D-1;i>=0;i--) 
	size[i] = L;
    stride[D-1] = 1;
    for(int i=D-2;i>=0;i--) 
	stride[i] = stride[i+1]*L;

    heff = new hType [N];

    mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a lattice (State) without setting the spin values
State::State(int dimension, int length, double disorder, int rngseed)  
{
  network=false;

    D = dimension;
    L = length;
    Z = 2*D;
    R = disorder;
    seed = rngseed;

    int i;
    N = 1;
    for(i=0;i<D;i++) 
	N *= L;     

    spin = new char [N];
    avalancheindex = new int [N];
    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    nUp = new int [N];

    Hext = 0;

    heff = new hType [N];

    nUp = new int [N];

    stride[D-1] = 1;
    for(int i=D-2;i>=0;i--) 
	stride[i] = stride[i+1]*L;

    for(int i=D-1;i>=0;i--) 
	size[i]=L;

    for(int i=0;i<Z;i++)
	neighborCoords[i] = new int[D];

    //for(int i=0; i<N; i++) 
    //spin[i] = DOWN; //spin[i] = UP;
    //M = -N;  
    //m = -1.0;
    //M = N;  
    //m = 1.0;
//    E0 = 0;
//    E  = 0;
    mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a lattice (State) with all spins = dir
State::State(int dimension, int length, double disorder, int rngseed, char dir)  
{
  network=false;

    D = dimension;
    L = length;
    Z = 2*D;
    R = disorder;
    seed = rngseed;

    int i;
    N = 1;
    for(i=0;i<D;i++) 
	N *= L;     

    spin = new char[N];
    avalancheindex = new int [N];
    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    nUp = new int [N];

    Hext = 0;
    heff = new hType[N];

    stride[D-1] = 1;
    for(int i=D-2;i>=0;i--) 
	stride[i] = stride[i+1]*L;

    for(int i=D-1;i>=0;i--) 
	size[i]=L;

    for(int i=0;i<Z;i++)
	neighborCoords[i] = new int[D];

    if(dir==1)
    {
      for(int i=0; i<N; i++) { 
	spin[i] = UP;
	nUp[i] = Z;
      }
      M = N;  
      m = 1.0;
    }
    else
      {
	for(int i=0; i<N; i++) {
	  spin[i] = DOWN;
	  nUp[i] = 0;
	}
 	M = -N;  
	m = -1.0;
    }

//    E0 = 0;
//    E  = 0;
    mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a lattice (State) with all spins values depending on external field H
State::State(int dimension, int length, double disorder, int rngseed, double H)  
{
  network=false;

    D = dimension;
    L = length;
    Z = 2*D;
    R = disorder;
    seed = rngseed;

    int i;
    N = 1;
    for(i=0;i<D;i++) 
	N *= L;     

    spin = new char [N];
    avalancheindex = new int [N];
    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    nUp = new int [N];

    Hext = H;
    heff = new hType [N];

    stride[D-1] = 1;
    for(int i=D-2;i>=0;i--) 
	stride[i] = stride[i+1]*L;

    for(int i=D-1;i>=0;i--) 
	size[i]=L;

    for(int i=0;i<Z;i++)
	neighborCoords[i] = new int[D];

    if(H>0)
    {
      for(int i=0; i<N; i++) {
	spin[i] = UP;
	nUp[i] = Z;
      }
	M = N;  
	m = 1.0;
    }
    else
    {
      for(int i=0; i<N; i++) { 
	spin[i] = DOWN;
	nUp[i] = 0;
      }
 	M = -N;  
	m = -1.0;
    }

    //E0 = 0;
    //E  = 0;

    mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a lattice (State) with all spins values depending on external field H
State::State(const State &X)  
{
  network=false;

    D = X.D;
    L = X.L;
    N = X.N;
    Z = 2*D;
    R = X.R;
    seed = X.seed;

    spin = new char [N];
    avalancheindex = new int [N];
    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    nUp = new int [N];

    Hext = X.Hext;
    heff = new hType [N];

    stride[D-1] = 1;
    for(int i=D-2;i>=0;i--) 
	stride[i] = stride[i+1]*L;

    for(int i=D-1;i>=0;i--) 
	size[i]=L;

    for(int i=0;i<Z;i++)
	neighborCoords[i] = new int[D];

    
    for(int i=0; i<N; i++) { 
	spin[i] = X.spin[i];
	nUp[i] = X.nUp[i];
    }

    M = X.M;
    m = X.m;
    E0 = X.E0;
    E  = X.E;

    mapping=false;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////// 
State::~State() 
{
  if(!network) {
    delete [] size; 
    delete [] stride;
    delete [] neighborLocs;
    for(int i=0;i<Z;i++)
	delete [] neighborCoords[i];
    delete [] neighborCoords;
  }

  delete [] spin; 
  delete [] avalancheindex;
  delete[] heff; 

  delete [] nUp;

}
////////////////////////////////////////////////////////////////////////////////////////////////////// 




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a network (State) with all spins = dir
State::State(int num, double disorder, int rngseed, char dir)  
{
  //cout << "test!!!"; //test 

  network = true; // this flag indicates that the system is a random network
  N = num;
  R = disorder;
  seed = rngseed;

  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int[N];
  Hext = 0;
  heff = new hType[N];

  if(dir==1) {
    for(int i=0; i<N; i++) {
      spin[i] = UP;
      nUp[i] = 0; // this should be K[i], to be calculated 
    }
    M = N;  
    m = 1.0;
  }
  else  {
    for(int i=0; i<N; i++) {
      spin[i] = DOWN;
      nUp[i] = 0; 
    }
    M = -N;  
    m = -1.0;
  }
    
  //    E0 = 0;
  //    E  = 0;
  mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a network (State) without spin values
State::State(int num, double disorder, int rngseed)  
{
  network = true; // this flag indicates that the system is a random network
  N = num;
  R = disorder;
  seed = rngseed;

  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int [N];

  Hext = 0;
  heff = new hType[N];

  //    E0 = 0;
  //    E  = 0;
  mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a real network (State) without spin values
State::State(int num, char dir)  
{
  network = true; // this flag indicates that the system is a random network. This is very important!! 

  N = num;
  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int [N];

  Hext = 0;
  heff = new hType[N];

  if(dir==1) {
    for(int i=0; i<N; i++) {
      spin[i] = UP;
      nUp[i] = 0; // should be K[i], to be calculated 
    }
    M = N;  
    m = 1.0;
  }
  else  {
    for(int i=0; i<N; i++) {
      spin[i] = DOWN;
      nUp[i] = 0; 
    }
    M = -N;  
    m = -1.0;
  }

  //    E0 = 0;
  //    E  = 0;
  mapping=false;



  // for single spin flip dynamics
  nUp = new int[N];
  


}
//////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a real network (State) without spin values
State::State(int num, double disorder, char dir)  
{
  network = true; // this flag indicates that the system is a random network. This is very important!! 

  N = num;
  R = disorder;
  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int [N];

  Hext = 0;
  heff = new hType[N];

  if(dir==1) {
    for(int i=0; i<N; i++) {
      spin[i] = UP;
      nUp[i] = 0; // should be K[i], to be calculated
    }
    M = N;  
    m = 1.0;
  }
  else  {
    for(int i=0; i<N; i++) {
      spin[i] = DOWN;
      nUp[i] = 0;
    }
    M = -N;  
    m = -1.0;
  }

  //    E0 = 0;
  //    E  = 0;
  mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a real network (State) without spin values
State::State(int num)  
{
  network = true; // this flag indicates that the system is a random network. This is very important!! 

  N = num;
  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int[N];

  Hext = 0;
  heff = new hType[N];

  //    E0 = 0;
  //    E  = 0;
  mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////////////////////////////// 
// construct a real network (State) without spin values
State::State(int num, double disorder)  
{
  network = true; // this flag indicates that the system is a random network. This is very important!! 

  N = num;
  R = disorder;
  spin = new char[N];
  avalancheindex = new int [N];
  nUp = new int [N];

  Hext = 0;
  heff = new hType[N];

  //    E0 = 0;
  //    E  = 0;
  mapping=false;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
