#include "State.h"

using namespace std;




void State::DFS()
{
  vector<bool> visited(N, false);
  
  AllComponents.clear();
  
  int maxsize  = 0; 
  int lccindex = 0; 
  int ncc  = -1;    
  for(int i=0; i<N; i++) {
    if(!visited[i]) {
      ncc++; 
      Nbl Component;
      DFS(i, visited, Component);
      AllComponents.push_back(Component);

      int size = Component.size();
      if(size>maxsize) {
	maxsize = size;
	lccindex = ncc;
      }

      
      
      
    }
  }
  

  Nlc = maxsize; 
  LC.clear();
  LC.resize(N);           


  Ncc = AllComponents.size();
  cout << "There are total " << Ncc << " connected components.\n";

  ccindex.clear();
  ccindex.resize(N);
  int cc = 0;	
  Nac = 0;
  for(vector<Nbl>::iterator q = AllComponents.begin(); q!=AllComponents.end(); q++) {
    for(Nbl_itr p = (*q).begin(); p!= (*q).end(); p++) {
      ccindex[*p]= cc;
    }
    int size = (*q).size();
    Nac += size;
    cc++;
  }
  Nac = (Nac-Nlc)/(Ncc-1);

  cc = 0;	
  int sumKlc = 0;    
  int Kofleader = 0; 
  for(vector<Nbl>::iterator q = AllComponents.begin(); q!=AllComponents.end(); q++) {
    if(cc==lccindex) { 
	
      int count=0;
      for(Nbl_itr p = (*q).begin(); p!= (*q).end(); p++, count++) {
	
	

	LClist.push_back(*p);    
	LC[(*p)] = 1;            

	if(K[(*p)] > Kofleader) {
	  leader = (*p);
	  Kofleader = K[leader];
	}
	  
	sumKlc += K[(*p)];
      }
      break;
    }
    else
      cc++;
  }
  Kmaxlc = Kofleader;
    
  
  Elc = sumKlc/2;
  mlc = Elc/(double)Nlc;
  cout << "\n Nlc = " << Nlc << " sumKlc= " << sumKlc << "  Elc=sumKlc/2= " << Elc 
       << " mlc = Elc/Nlc= " << mlc << endl;

  
  
  
  


   


}


 
void State::DFS(int u, vector<bool>& visited, Nbl&  Component)
{
  visited[u] = true;
  Component.push_back(u);
    
  
  
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    
    int v = (*p);
    
    
    if(!visited[v])
      DFS(v, visited, Component);
  }
}
