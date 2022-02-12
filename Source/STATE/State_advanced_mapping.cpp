#include "State.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////// 
// map the random field to a network file (DIMACS format). This function demonstrates the mapping in 
// details though it is not used in practice (We do the mapping directly in parser.cpp.). This 
// following function is based on Mapping(), it will map the RFIM to a network with some fixed spin 
// values. This is a step of using the earlier solution to speed up the program. 
//////////////////////////////////////////////////////////////////////////////////////////////////////
//DIR=1, freeze UP spins; DIR=-1, freeze DOWN spins; 

void State::Mapping_Frozenspins(char FZdir) 
{
  static int count=0;

  ofstream outfile1("file1", ios_base::out);   

  // set up the map between the free spin's location and the node index
  Node_Map_Loc(FZdir);

  int nEdge = 0;

  if(!network) { // if for a hypercubic lattice system 
    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {
	    GetNeighbors(i);

	    int nfjbig    = 0;   // count how many free neighbors with location index bigger than i
	    int nfjsmall  = 0;   // count how many free neighbors with location index smaller than i
	    int nFZnb = 0;       // the number of spin i's frozen neighbors
	    for(int j=0; j<Z; j++)
	      {
		if(spin[neighborLocs[j]] != FZdir) 
		  {
		    if(neighborLocs[j] > i) 
		      {
			nfjbig++;
			outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << NodeMap[neighborLocs[j]] << ' ' << CC << endl;
			nEdge++;
		      }
		    else
		      nfjsmall++;
		  }
		else
		  nFZnb++;		    
	      }

	    //hType Wi = -2*(heff[i] + (2*FZdir-1)*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field, // FZdir = 1 or 0
	    hType Wi = -2*(heff[i] + FZdir*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field,         // FZdir = 1 or -1

	    //*************************     TAKE CARE OF THE FOLLOWING COMMENTS         *******************************
            // here we have to consider only those bonds which connect i to FREE neighbor spins (not a frozen spin)
	    // Wi = -2 (hi +/- J sum_nnFZj ) - 1/2 * sum_freej (Cij-Cji) 
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * [sum_(freej>i) Cij - sum_(freej<i) Cji]
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * CC * [sum_(freej>i)  - sum_(freej<i) ]           here CC=4J
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * CC * [nfjbig - nfjsmall]
	    //    = -2 (hi +/- J sum_nnFZj ) - 2J * [nfjbig - nfjsmall]
	    // Here, + for FZdir = UP   = 1 ---> (2*FZdir-1) = +1
	    //       - for FZdir = DOWN = 0 ---> (2*FZdir-1) = -1
            //*********************************************************************************************************

	    if(Wi>0)
	      {
		outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << Nfs+1 << ' ' << Wi << endl;
		nEdge++;
	      }
	    else //if (Wi<=0)  // including the 0-capacity arc to avoid strange segmentation fault sometimes. Y.L 06/15/06
	      {
		outfile1<< 'a' << ' ' << 0  << ' ' << NodeMap[i] << ' ' << -Wi << endl;
		nEdge++;
	      }
    
	  }
      }// end of looping over node 
  }// end of if a lattice system
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  else { // else for a network system 
    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {

	    int nfjbig    = 0;   // count how many free neighbors with location index bigger than i
	    int nfjsmall  = 0;   // count how many free neighbors with location index smaller than i
	    int nFZnb = 0;       // the number of spin i's frozen neighbors

	    hType Wi = -2 * heff[i];
	    for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) 
	      {
		int j = *p;
		if(spin[j] != FZdir) 
		  {
		    if(j > i) 
		      {
			nfjbig++;
			hType Cij = CC * edgeweight(i, j);
			outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << NodeMap[j] << ' ' << Cij << endl;
			nEdge++;
			Wi -= 0.5*Cij;
		      }
		    else
		      {
			nfjsmall++;
			hType Cji = CC * edgeweight(j, i);
			Wi += 0.5*Cji;
		      }
		  }
		else {
		  nFZnb++;
		  hType Jij = J * edgeweight(i, j); 
		  //Wi -= 2 * (2*FZdir-1) * Jij;// FZdir = 1 or 0
		  Wi -= 2 * FZdir * Jij;        // FZdir = 1 or -1
		}		    
	      }


	    //hType Wi = -2*(heff[i] + (2*FZdir-1)*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field, 
	    //*************************     TAKE CARE OF THE FOLLOWING COMMENTS         *******************************
            // here we have to consider only those bonds which connect i to FREE neighbor spins (not a frozen spin)
	    // Wi = -2 (hi +/-  sum_nnFZj Jij) - 1/2 * sum_freej (Cij-Cji) 
	    //    = -2 (hi +/-  sum_nnFZj Jij) - 1/2 * [sum_(freej>i) Cij - sum_(freej<i) Cji]
	    // Here, + for FZdir = UP   = 1 ---> (2*FZdir-1) = +1
	    //       - for FZdir = DOWN = 0 ---> (2*FZdir-1) = -1
            //*********************************************************************************************************

	    if(Wi>0)
	      {
		outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << Nfs+1 << ' ' << Wi << endl;
		nEdge++;
	      }
	    else //if (Wi<=0)  // including the 0-capacity arc to avoid strange segmentation fault sometimes. Y.L 06/15/06
	      {
		outfile1<< 'a' << ' ' << 0  << ' ' << NodeMap[i] << ' ' << -Wi << endl;
		nEdge++;
	      }
    
	  }
      }// end of looping over node 
  }// end of if a network system





  outfile1.close();

  ofstream outfile2("file2", ios_base::out);
  outfile2<< "p  max " << Nfs + 2 << ' ' << nEdge <<endl; // pay attention to the # of nodes
  outfile2<< "n  0  s" <<endl;
  outfile2<< "n  " << Nfs + 1 << "  t" <<endl;                   
  outfile2.close();


  char cmd[256];
  sprintf(cmd,"cat file2 file1 > graph-%d.inp",count++);
  system(cmd);
  system("rm file*");
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
/* About the mapping, for example D=1, N=14, there are nUP=9 frozen UP spins:
   1D spin chain:        1 1 0 1 1 1 0 0 1 1  1  1  0  0
   spin array index:     0 1 2 3 4 5 6 7 8 9  10 11 12 13  (spin location)
   node index (old):     1 2 3 4 5 6 7 8 9 10 11 12 13 14  (just loc + 1)
   node index (new):         1       2 3            4  5

   therefore, there is a one-one mapping between the new node index and the
   spin array index (location) as following:

   loc <--->  node <---> loc
   2          1          2
   6          2          6   
   7          3          7
   12         4          12
   13         5          13 
   .          .          .
   .          .          .
   .          .          .
   ?         Nfs         ?

   so source s= 0, sink t= Nfs+1, n=Nfs+2, m=?

   Note that in the old code (without considering the frozen spin), the
   mapping between the spin location and the node index is simply:
   loc <---> node
   0          1
   .          .    
   .          .    
   .          .    
   i         i+1
   .          .    
   .          .    
   .          .    
   N-1         N

   so source s= 0, sink t= N+1, n=N+2, m=?
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::Node_Map_Loc(char FZdir)
{
  int tmp;
  if(!mapping) // set up the map between the free spin's location and the node index
    {
      // get the number of free sites. if we freeze UP spins, then Nfs=nDOWN=(N-M)/2; otherwise Nfs=nUP=(N+M)/2;
      if (FZdir==1) Nfs = (N-M)/2; 
      else 	    Nfs = (N+M)/2 ; 

      //cout << "I am here with N= " << N << "; M= " << M << "; Nfs= " << Nfs << endl; //debug
      // Note that if M has not been calculated, then we may have trouble, e.g., N= 20611; M= 1625051553; Nfs= -812515471
      // This causes error: libc++abi.dylib: terminating with uncaught exception of type std::length_error: vector
      LocVec.resize(1+Nfs); 

      tmp = 1;
      for(int i=0; i<N; i++)
	{ 
	  if(spin[i] != FZdir)  
	    {
	      NodeMap[i]  = tmp;
	      LocVec[tmp] = i;
	      tmp++;
	    }
	}
      mapping=true;
    }

  

}
//////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::Mapping_Frozenspins(char FZdir,long &n, long &m, long &s, long &t)
{
  //cout << "\n Mapping with frozen spins " << endl; //debug

  // set up the map between the free spin's location and the node index
  Node_Map_Loc(FZdir);

    
  int nEdge = 0;


  if(!network) { // if for a hypercubic lattice system 
 
    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir)
	  {
	    GetNeighbors(i);

	    int nfjbig    = 0;   // count how many free neighbors with location index bigger than i
	    int nfjsmall  = 0;   // count how many free neighbors with location index smaller than i
	    int nFZnb = 0;       // the number of spin i's frozen neighbors
	    for(int j=0; j<Z; j++)
	      {
		if(spin[neighborLocs[j]] != FZdir) 
		  {
		    if(neighborLocs[j] > i) 
		      {
			nfjbig++;
			myarc arc(NodeMap[i], NodeMap[neighborLocs[j]], CC);
			ArcList.push_back(arc);
			nEdge++;
			sum_Cij += CC;
		      }
		    else
		      nfjsmall++;
		  }
		else
		  nFZnb++;		    
	      }

	    //hType Wi = -2*(heff[i] + (2*FZdir-1)*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field,// FZdir = 1 or 0
	    hType Wi = -2*(heff[i] + FZdir*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field,        // FZdir = 1 or -1
 
	    //*************************     TAKE CARE OF THE FOLLOWING COMMENTS         *******************************
	    // here we have to consider only those bonds which connect i to FREE neighbor spins (not a frozen spin)
	    // Wi = -2 (hi +/- J sum_nnFZj ) - 1/2 * sum_freej (Cij-Cji) 
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * [sum_(freej>i) Cij - sum_(freej<i) Cji]
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * CC * [sum_(freej>i)  - sum_(freej<i) ]           here CC=4J
	    //    = -2 (hi +/- J sum_nnFZj ) - 1/2 * CC * [nfjbig - nfjsmall]
	    //    = -2 (hi +/- J sum_nnFZj ) - 2J * [nfjbig - nfjsmall]
	    // Here, + for FZdir = UP   = 1 ---> (2*FZdir-1) = +1
	    //       - for FZdir = DOWN = 0 ---> (2*FZdir-1) = -1
	    //*********************************************************************************************************

	    if(Wi>0)
	      {
		myarc arc(NodeMap[i], Nfs+1, Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_pos++;
		sum_absWi += Wi;

	      }
	    else //if (Wi<=0)  // including the 0-capacity arc to avoid strange segmentation fault sometimes. Y.L 06/15/06
	      {
		myarc arc(0,NodeMap[i], -Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_neg++;
		sum_absWi -= Wi;

	      }
    
	  }

	//cout << "  # of positive Wi = " << num_Wi_pos << "; # of negative Wi = " << num_Wi_neg << ' '; 
	// Note that if all Wi's are positive, then the source node will be isolated.
	// Similarly,if all Wi's are negative, then the sink   node will be isolated.   
	// This will cause error: "line 0 of input - source or sink doesn't have incident arcs."
	// To avoid this issue, we need to consider an edge connecting source and sink directly:
	
	hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
	//cout << Cst << endl; //test
	//myarc arc(0, N+1, Cst);
	myarc arc(0, Nfs+1, Cst);
	ArcList.push_back(arc);
	nEdge++;
	

      }
  }// end of a hypercubic lattice system
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  else { // else for a network system 

   int num_Wi_pos = 0;
    int num_Wi_neg = 0;
 
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {

	    int nfjbig   = 0;   // count how many free neighbors with location index bigger than i
	    int nfjsmall = 0;   // count how many free neighbors with location index smaller than i
	    int nFZnb = 0;      // the number of spin i's frozen neighbors

	    hType Wi = -2 * heff[i];
	    for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) 
	      {
		int j = *p;
		if(spin[j] != FZdir) 
		  {
		    if(j > i) 
		      {
			nfjbig++;
			hType Cij = CC * edgeweight(i, j);
			myarc arc(NodeMap[i], NodeMap[j], Cij);
			ArcList.push_back(arc);
			nEdge++;
			Wi -= 0.5*Cij;
			sum_Cij += Cij;

		      }
		    else
		      {
			nfjsmall++;
			hType Cji = CC * edgeweight(j, i);
			Wi += 0.5*Cji;
		      }
		  }
		else {
		  nFZnb++;
		  hType Jij = J * edgeweight(i, j); 
		  //Wi -= 2 * (2*FZdir-1) * Jij; // FZdir = 1 or 0
		  Wi -= 2 * FZdir * Jij;         // FZdir = 1 or -1  
		}		    

	      }

	    //hType Wi = -2*(heff[i] + (2*FZdir-1)*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  // an auxiliary field, 
	    //*************************     TAKE CARE OF THE FOLLOWING COMMENTS         *******************************
            // here we have to consider only those bonds which connect i to FREE neighbor spins (not a frozen spin)
	    // Wi = -2 (hi +/-  sum_nnFZj Jij) - 1/2 * sum_freej (Cij-Cji) 
	    //    = -2 (hi +/-  sum_nnFZj Jij) - 1/2 * [sum_(freej>i) Cij - sum_(freej<i) Cji]
	    // Here, + for FZdir = UP   = 1 ---> (2*FZdir-1) = +1
	    //       - for FZdir = DOWN = 0 ---> (2*FZdir-1) = -1
            //*********************************************************************************************************

	    if(Wi>0)
	      {
		myarc arc(NodeMap[i], Nfs+1, Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_pos++;
		sum_absWi += Wi;
	      }
	    else //if (Wi<=0)  // including the 0-capacity arc to avoid strange segmentation fault sometimes. Y.L 06/15/06
	      {
		myarc arc(0,NodeMap[i], -Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_neg++;
		sum_absWi -= Wi;
	      }
    
	  }

	//cout << "  # of positive Wi = " << num_Wi_pos << "; # of negative Wi = " << num_Wi_neg << ' '; 
	hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
	//cout << ' ' << Cst << endl; //test 
	//myarc arc(0, N+1, Cst);
	myarc arc(0, Nfs+1, Cst);
	ArcList.push_back(arc);
	nEdge++;
   
      }// end of looping over node 
  }// end of if a network system




  n = Nfs + 2; 
  m = nEdge; 
  s = 0;
  t = Nfs + 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
