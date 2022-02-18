#include "State.h"

using namespace std;


void State::Mapping()
{
  static int count=0;

  ofstream outfile1("file1", ios_base::out);

  int nEdge = 0;

  if(network==false) { 
    
    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      GetNeighbors(i);
      
      int njbig=0; 
      for(int j=0; j<Z; j++) {
	if(neighborLocs[j]>i) {
	  njbig++;
	  outfile1 << 'a' << ' '<< 1+i << ' ' << 1+neighborLocs[j] << ' ' << CC << endl;
	  nEdge++;
	  sum_Cij += CC;
	}
      }
      
      hType Wi = -2*heff[i] - CC*(njbig-D);  
      
      if(Wi>0){
	outfile1 << 'a' << ' ' << 1+i << ' ' << N+1 << ' ' << Wi << endl;
	nEdge++;
	num_Wi_pos++;
	sum_absWi += Wi;
      }
      else {
	outfile1 << 'a' << ' ' <<  0  << ' ' << 1+i << ' ' << -Wi << endl;
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }

    cout << " # of positive Wi = " << num_Wi_pos 
    	 << ";# of negative Wi = " << num_Wi_neg << endl; 
    
    
    
    
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
  
    outfile1 << 'a' << ' ' << 0 << ' ' << N+1 << ' ' << Cst << endl;
    nEdge++;
   
  } 
  

  
  else { 

    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      hType Wi = -2*heff[i];
      
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
      else {
	outfile1 << 'a' << ' ' <<  0  << ' ' << 1+i << ' ' << -Wi << endl;
	nEdge++;
      }
    }

    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
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






void State::Mapping(long &n, long &m, long &s, long &t)
{
 
  int nEdge = 0;
 
  if(network==false) { 

    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      GetNeighbors(i);
      int njbig=0; 
      
      for(int j=0; j<Z; j++) 	  {
	if(neighborLocs[j]>i)    {
	  njbig++;
	  myarc arc(1+i, 1+neighborLocs[j], CC);
	  ArcList.push_back(arc);
	  nEdge++;
	  sum_Cij += CC;
	}
      }
      
      hType Wi = -2*heff[i] - CC*(njbig-D);  
    
      if(Wi>0){
	myarc arc(1+i, N+1, Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_pos++;
	sum_absWi += Wi;
      }
      else {
	myarc arc(0, 1+i, -Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }
    
    
    
    
    
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
    
    myarc arc(0, N+1, Cst);
    ArcList.push_back(arc);
    nEdge++;
    

  } 
  


  
  else { 
   int num_Wi_pos = 0;
    int num_Wi_neg = 0;
 
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) {
      hType Wi = -2*heff[i];
      
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = *p;
	if(i<j) {
	  

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
      else {
	myarc arc(0, 1+i, -Wi);
	ArcList.push_back(arc);
	nEdge++;
	num_Wi_neg++;
	sum_absWi -= Wi;
      }
    }

    
    hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
    
    myarc arc(0, N+1, Cst);
    ArcList.push_back(arc);
    nEdge++;
    
  } 
 


  n = N + 2;
  m = nEdge;
  s = 0;
  t = N + 1;
}


