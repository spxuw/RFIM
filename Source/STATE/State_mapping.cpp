#include "State.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////// 
// map the random field to a network file (DIMACS format). This function demonstrates the mapping in details 
// though it is not used in practice (We do the mapping directly in parser.cpp.).

/* Attention: finally in the DIMACS file,  
   node   0   represents the source s
   node   1   represents the 1st node (with spin location index 0)
   node   2   represents the 2nd node (with spin location index 1)
   .............
   node   N   represents the Nth node (with spin location index N-1)
   node   N+1 represents the sink   t 

   Also note that for a special case: D=2,L=2 with PBC, the nearest nerighbors are double-counted.
   0 1 
   2 3 
   Considering PBC, spin 0 has nn: 1,2,1,2. spin 1 has nn: 0,3,0,3. and so on. Consequently, there will 
   affect the arc descriptions. */ 

/* About the calculation of the auxiliary field. 

   Note that: Cij = 0 if i>=j; 4Jij if i<j.   

   Auxiliary field (see my thesis Eq.A.14)
   Wi = -2*heff[i] - 1/2 * sum_j [Cij - Cji]
   = -2*heff[i] - 1/2 * [sum_(j>i) Cij - sum_(j<i) Cji]

   For lattice: Z=2D, Jij=1, Cij=CC=4J
   Wi = -2*heff[i] - 1/2 * CC * [sum_(j>i)  - sum_(j<i) ]     here CC=4J
   = -2*heff[i] - 1/2 * CC * [njbig - (Z-njbig)]           here Z=2D  
   = -2*heff[i] - CC*(njbig-D)

   For network, Cij = CC * Jij  
   Wi = -2*heff[i] - 1/2 * [sum_(j>i) Cij - sum_(j<i) Cji]

*/



//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::Mapping()
{
  static int count=0;

  ofstream outfile1("file1", ios_base::out);

  int nEdge = 0;

  if(network==false) { // if we have a hypercubic lattice system
    
    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      GetNeighbors(i);
      
      int njbig=0; //count how many neighbors with location index bigger than i
      for(int j=0; j<Z; j++) {
	if(neighborLocs[j]>i) {
	  njbig++;
	  outfile1 << 'a' << ' '<< 1+i << ' ' << 1+neighborLocs[j] << ' ' << CC << endl;
	  nEdge++;
	  sum_Cij += CC;
	}
      }
      
      hType Wi = -2*heff[i] - CC*(njbig-D);  // an auxiliary field, 
      
      if(Wi>0){
	outfile1 << 'a' << ' ' << 1+i << ' ' << N+1 << ' ' << Wi << endl;
	nEdge++;
	num_Wi_pos++;
	sum_absWi += Wi;
      }
      else {// if (Wi<0)
	outfile1 << 'a' << ' ' <<  0  << ' ' << 1+i << ' ' << -Wi << endl;
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }

    cout << " # of positive Wi = " << num_Wi_pos 
    	 << ";# of negative Wi = " << num_Wi_neg << endl; 
    // Note that if all Wi's are positive, then the source node will be isolated.
    // Similarly,if all Wi's are negative, then the sink   node will be isolated.   
    // This will cause error: "line 0 of input - source or sink doesn't have incident arcs."
    // To avoid this issue, we need to consider an edge connecting source and sink directly:
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
  
    outfile1 << 'a' << ' ' << 0 << ' ' << N+1 << ' ' << Cst << endl;
    nEdge++;
   
  } // end of if a lattice
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  else { // for a network system 

    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      hType Wi = -2*heff[i];// an auxiliary field
      
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = *p;
	if(i<j) {
	  hType Cij = CC * edgeweight(i, j);
	  outfile1 << 'a' << ' '<< 1+i << ' ' << 1+j << ' ' << Cij << endl;
	  nEdge++;
	  Wi -= 0.5*Cij;
	  sum_Cij += Cij;
	}
	else if (j<i) {  
	  hType Cji = CC * edgeweight(j, i);  
	  Wi += 0.5*Cji;  
	}
      }
      
      if(Wi>0){
	outfile1 << 'a' << ' ' << 1+i << ' ' << N+1 << ' ' << Wi << endl;
	nEdge++;
      }
      else {// if (Wi<0)
	outfile1 << 'a' << ' ' <<  0  << ' ' << 1+i << ' ' << -Wi << endl;
	nEdge++;
      }
    }

    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
    outfile1 << 'a' << ' ' << 0 << ' ' << N+1 << ' ' << Cst << endl;
     nEdge++;

  }



  outfile1.close();

  ofstream outfile2("file2", ios_base::out);
  outfile2<< "p  max " << N+2 << ' ' << nEdge <<endl;
  outfile2<< "n  0  s" <<endl;
  outfile2<< "n  " << N+1 << "  t" <<endl;                   
  outfile2.close();


  char cmd[256];
  sprintf(cmd,"cat file2 file1 > graph-true-%d.inp",count++);
  system(cmd);
  system("rm file*");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::Mapping(long &n, long &m, long &s, long &t)
{
 
  int nEdge = 0;
 
  if(network==false) { // if we have a hypercubic lattice system

    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      GetNeighbors(i);
      int njbig=0; //count how many neighbors with location index bigger than i
      
      for(int j=0; j<Z; j++) 	  {
	if(neighborLocs[j]>i)    {
	  njbig++;
	  myarc arc(1+i, 1+neighborLocs[j], CC);
	  ArcList.push_back(arc);
	  nEdge++;
	  sum_Cij += CC;
	}
      }
      
      hType Wi = -2*heff[i] - CC*(njbig-D);  // an auxiliary field, 
    
      if(Wi>0){
	myarc arc(1+i, N+1, Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_pos++;
	sum_absWi += Wi;
      }
      else {// if (Wi<0)
	myarc arc(0, 1+i, -Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }
    //cout << "  # of positive Wi = " << num_Wi_pos << "; # of negative Wi = " << num_Wi_neg << endl; 
    // Note that if all Wi's are positive, then the source node will be isolated.
    // Similarly,if all Wi's are negative, then the sink   node will be isolated.   
    // This will cause error: "line 0 of input - source or sink doesn't have incident arcs."
    // To avoid this issue, we need to consider an edge connecting source and sink directly:
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
    //cout << Cst << endl;
    myarc arc(0, N+1, Cst);
    ArcList.push_back(arc);
    nEdge++;
    

  } // end of if a lattice
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  else { // for a network system 
   int num_Wi_pos = 0;
    int num_Wi_neg = 0;
 
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      hType Wi = -2*heff[i];// an auxiliary field
      
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = *p;
	if(i<j) {
	  //cout <<  edgeweight(i, j) << ','; //test 

	  hType Cij = CC * edgeweight(i, j); 
	  myarc arc(1+i, 1+j, Cij);
	  ArcList.push_back(arc);
	  nEdge++;
	  Wi -= 0.5*Cij;
	  sum_Cij += Cij;
	}
	else if (j<i) {  
	  hType Cji = CC * edgeweight(j, i);  
	  Wi += 0.5*Cji;  
	}
      }
      
      if(Wi>0){
	myarc arc(1+i, N+1, Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_pos++;
	sum_absWi += Wi;
      }
      else {// if (Wi<0)
	myarc arc(0, 1+i, -Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }

    //cout << "  # of positive Wi = " << num_Wi_pos << "; # of negative Wi = " << num_Wi_neg << endl; 
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; // See Eq.(A.16) in my thesis
    //cout << Cst << endl; 
    myarc arc(0, N+1, Cst);
    ArcList.push_back(arc);
    nEdge++;
    
  } // end of if a network 
 


  n = N + 2;
  m = nEdge;
  s = 0;
  t = N + 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////

