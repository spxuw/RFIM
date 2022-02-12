#include "State.h"

using namespace std;

// the following DFS is designed to find the largest connected component in a directed or undirected graph
// for directed graph, the DFS is performed by searching the A_list, rather than Aout_list, of each node.
// So it cannot give the topological order of nodes.
void State::DFS()
{
  vector<bool> visited(N, false);
  //vector<Nbl>  AllComponents; // all the connected components
  AllComponents.clear();
  
  int maxsize  = 0; // size of the largest connected component
  int lccindex = 0; // index of the largest connected component
  int ncc  = -1;    // count the number of connected component
  for(int i=0; i<N; i++) {
    if(!visited[i]) {
      ncc++; // find another connected component
      Nbl Component;
      DFS(i, visited, Component);
      AllComponents.push_back(Component);

      int size = Component.size();
      if(size>maxsize) {
	maxsize = size;
	lccindex = ncc;
      }

      //cout << "\n Component with size ("<< size  << "): "; //test
      //for(Nbl_itr p = Component.begin(); p!= Component.end(); p++) 
      //cout << (*p) << ' ';
    }
  }
  //cout << "\nThe largest connected component (label " << lccindex << ") has size = " << maxsize; //test

  Nlc = maxsize; // return the size of the largest connected component
  LC.clear();
  LC.resize(N);           // this vector labeling the node whether it belongs to the largest connected component


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
  int sumKlc = 0;    // sum of k for the largest connected component
  int Kofleader = 0; // degree of the leader (the node with the highest degree in LCC)
  for(vector<Nbl>::iterator q = AllComponents.begin(); q!=AllComponents.end(); q++) {
    if(cc==lccindex) { // if this is the largest component
	
      int count=0;
      for(Nbl_itr p = (*q).begin(); p!= (*q).end(); p++, count++) {
	//cout << (*p) << ' ' << K[(*p)] << endl; //debug
	//cout << (*p) << ' '; //debug

	LClist.push_back(*p);    // push it into the LClist
	LC[(*p)] = 1;            // label it 

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
    
  // calculate the number of edges in the largest connected component
  Elc = sumKlc/2;
  mlc = Elc/(double)Nlc;
  cout << "\n Nlc = " << Nlc << " sumKlc= " << sumKlc << "  Elc=sumKlc/2= " << Elc 
       << " mlc = Elc/Nlc= " << mlc << endl;

  //cout << "\nLC[] = ";
  //for(int i=0; i<N; i++) 
  //cout << LC[i] << ' '; //test
  //cout << " with leader node " << leader << ".\n";   //test


   


}


 
void State::DFS(int u, vector<bool>& visited, Nbl&  Component)
{
  visited[u] = true;
  Component.push_back(u);
    
  //cout << "u = " << u ; //test
  // loop over node u's adjacent nodes
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    //for(Nbl_itr p = pseudoA[u].begin(); p!= pseudoA[u].end(); p++) { // test for directed graph
    int v = (*p);
    //cout << "v= " << v << endl ; //test
    //printf(" v= %d \n", v); //test
    if(!visited[v])
      DFS(v, visited, Component);
  }
}
