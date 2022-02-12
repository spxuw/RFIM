#include "State.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////// 
// initialize the lattice (State) with all spins = dir
void State::Initialize_withdir(char dir)  
{
    if(dir==1)
    {
	for(int i=0; i<N; i++) 
	    spin[i] = UP;
	M = N;  
	m = 1.0;
    }
    else
    {
	for(int i=0; i<N; i++) 
	    spin[i] = DOWN;
 	M = -N;  
	m = -1.0;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////// 
// initialize the lattice (State) with spins UP if H>0 
void State::Initialize_fromfield(double H)  
{
    if(H<=0.0)
    {
	for(int i=0; i<N; i++) 
	    spin[i] = DOWN;
 	M = -N;  
	m = -1.0;
    }
    else
    {
	for(int i=0; i<N; i++) 
	    spin[i] = UP;
	M = N;  
	m = 1.0;
    }

}
//////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////// 
// initialize State from the GS file: ./data/GS-jobID/k
void State::Initialize_fromfile(int k)  
{
    char fname[256]; 
    if(!network)
      sprintf(fname,"./data/GS-D%d-L%d-R%.3lf-seed%d/%d",D,L,R,seed,k);
    else
      sprintf(fname,"./data/GS-N%d-c%e-R%.3lf-seed%d/%d",N,c,R,seed,k);


    ifstream infile(fname, ios::in|ios::binary);
    if (!infile)  
    {
	cout << "cannot open file:  " << fname << "for read\n";
	exit(2);
    }
    
    //infile >> M;          // M must be known to calculate the number of free sites (Nfs)
    infile.read(reinterpret_cast<char *>(&M),sizeof(int));

    m = M/(double)N;

    //for(int i=0;i<N;i++)
    //infile >> spin[i];
    infile.read(spin,N*sizeof(char));
    
    
    infile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////////////////////////////////// 
// initialize State from the GS file: ./data/GS-jobID/k
void State::Initialize_fromfile(string file, int k)  
{
    char fname[256]; 
    sprintf(fname,"./data/GS-%s/%d",file.c_str(),k);
    
    ifstream infile(fname, ios::in|ios::binary);
    if (!infile)  
    {
	cout << "cannot open file:  " << fname << "for read\n";
	exit(2);
    }
    
    //infile >> M;          // M must be known to calculate the number of free sites (Nfs)
    infile.read(reinterpret_cast<char *>(&M),sizeof(int));

    //cout << "\nM= " << M << endl; //debug 

    m = M/(double)N;

    //for(int i=0;i<N;i++)
    //infile >> spin[i];
    infile.read(spin,N*sizeof(char));
    
    
    infile.close();
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////// 
// Reset for the calculation of ground state at different external field
void State::Reset()
{
    mapping = false; 
    ArcList.clear();
}
////////////////////////////////////////////////////////////////////////////////////////////////////// 
