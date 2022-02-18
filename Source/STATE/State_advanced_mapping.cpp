#include "State.h"

using namespace std;

void State::Mapping_Frozenspins(char FZdir) 
{
  static int count=0;

  ofstream outfile1("file1", ios_base::out);   

  
  Node_Map_Loc(FZdir);

  int nEdge = 0;

  if(!network) { 
    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {
	    GetNeighbors(i);

	    int nfjbig    = 0;   
	    int nfjsmall  = 0;   
	    int nFZnb = 0;       
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

	    
	    hType Wi = -2*(heff[i] + FZdir*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  

	    
            
	    
	    
	    
	    
	    
	    
	    
            

	    if(Wi>0)
	      {
		outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << Nfs+1 << ' ' << Wi << endl;
		nEdge++;
	      }
	    else 
	      {
		outfile1<< 'a' << ' ' << 0  << ' ' << NodeMap[i] << ' ' << -Wi << endl;
		nEdge++;
	      }
    
	  }
      }
  }
  



  
  else { 
    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {

	    int nfjbig    = 0;   
	    int nfjsmall  = 0;   
	    int nFZnb = 0;       

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
		  
		  Wi -= 2 * FZdir * Jij;        
		}		    
	      }


	    
	    
            
	    
	    
	    
	    
            

	    if(Wi>0)
	      {
		outfile1<< 'a' << ' ' << NodeMap[i] << ' ' << Nfs+1 << ' ' << Wi << endl;
		nEdge++;
	      }
	    else 
	      {
		outfile1<< 'a' << ' ' << 0  << ' ' << NodeMap[i] << ' ' << -Wi << endl;
		nEdge++;
	      }
    
	  }
      }
  }





  outfile1.close();

  ofstream outfile2("file2", ios_base::out);
  outfile2<< "p  max " << Nfs + 2 << ' ' << nEdge <<endl; 
  outfile2<< "n  0  s" <<endl;
  outfile2<< "n  " << Nfs + 1 << "  t" <<endl;                   
  outfile2.close();


  char cmd[256];
  sprintf(cmd,"cat file2 file1 > graph-%d.inp",count++);
  system(cmd);
  system("rm file*");
    
}







void State::Node_Map_Loc(char FZdir)
{
  int tmp;
  if(!mapping) 
    {
      
      if (FZdir==1) Nfs = (N-M)/2; 
      else 	    Nfs = (N+M)/2 ; 

      
      
      
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






void State::Mapping_Frozenspins(char FZdir,long &n, long &m, long &s, long &t)
{
  

  
  Node_Map_Loc(FZdir);

    
  int nEdge = 0;


  if(!network) { 
 
    int num_Wi_pos = 0;
    int num_Wi_neg = 0;
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir)
	  {
	    GetNeighbors(i);

	    int nfjbig    = 0;   
	    int nfjsmall  = 0;   
	    int nFZnb = 0;       
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

	    
	    hType Wi = -2*(heff[i] + FZdir*J*nFZnb) - 2*J*(nfjbig - nfjsmall);  
 
	    
	    
	    
	    
	    
	    
	    
	    
	    
	    

	    if(Wi>0)
	      {
		myarc arc(NodeMap[i], Nfs+1, Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_pos++;
		sum_absWi += Wi;

	      }
	    else 
	      {
		myarc arc(0,NodeMap[i], -Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_neg++;
		sum_absWi -= Wi;

	      }
    
	  }

	
	
	
	
	
	
	hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
	
	
	myarc arc(0, Nfs+1, Cst);
	ArcList.push_back(arc);
	nEdge++;
	

      }
  }
  


  
  else { 

   int num_Wi_pos = 0;
    int num_Wi_neg = 0;
 
    hType sum_Cij = 0;
    hType sum_absWi = 0;

    for(int i=0; i<N; i++) 
      {
	if(spin[i] != FZdir) 
	  {

	    int nfjbig   = 0;   
	    int nfjsmall = 0;   
	    int nFZnb = 0;      

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
		  
		  Wi -= 2 * FZdir * Jij;         
		}		    

	      }

	    
	    
            
	    
	    
	    
	    
            

	    if(Wi>0)
	      {
		myarc arc(NodeMap[i], Nfs+1, Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_pos++;
		sum_absWi += Wi;
	      }
	    else 
	      {
		myarc arc(0,NodeMap[i], -Wi);
		ArcList.push_back(arc);
		nEdge++;
		num_Wi_neg++;
		sum_absWi -= Wi;
	      }
    
	  }

	
	hType Cst = -0.25*sum_Cij - 0.5*sum_absWi; 
	
	
	myarc arc(0, Nfs+1, Cst);
	ArcList.push_back(arc);
	nEdge++;
   
      }
  }




  n = Nfs + 2; 
  m = nEdge; 
  s = 0;
  t = Nfs + 1;
}

