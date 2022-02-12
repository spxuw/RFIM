#include "State.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////////
// i--j1--k
// i--j2--k
// put the common neighbors of node i and node k into a list CN
int State::GetHighestscoringCommonNeighbors(int i, int k)
{
  double Smax = -N;
  int j0 = 0;

  for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
    int j = (*p); 
    for(Nbl_itr q = A[j].begin(); q!= A[j].end(); q++) {
      if((*q)==k) {
	if(h[j]>Smax) {
	  Smax = h[j];
	  j0 = j;
	}
      }
    }
  }
  return j0;
}
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
/* Highscoring Module Search using Peiling Jia's greedy algorithm: each node is a seed gene  
   Note that if the next-nearest neighbor k of node i is included, 
   then we also need to decide which nearest neighbor that lead us to k should be included
   (because they could be multiple nodes connecting node i with node k)
   Here, we just choose the one with the highest score.

   i--j1--k
   i--j2--k
   If node j1 and j2 have different weights, edges (i,j1), (j1,k), (i,j2), (j2,k) also have different weights.
   What shall we do? Which path shall we choose?



*/
void State::GreedyGrowth(Nbl& NL, double r)
{
  double S0, S, S0new, Snew;
  CalModuleScore(NL, S0, S);

  vector<bool> visited(N, false);
  for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) {
    int i = (*p);
    if(!visited[i]) {
      
      int size0 = NL.size();
      //check its nearest neighbor 
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = (*p); 
	
	if(!find(NL, j)) {
	  NL.push_back(j);
	  
	  CalModuleScore(NL, S0new, Snew); 
	  
	  if(Snew>S*(1+r))
	    S = Snew;
	  //if(S0new>S0*(1+r)) 
	  //S0 = S0new;
	  else  
	    NL.remove(j); 
	}
      }// end of checking nearest neighbor of node i
      
      // if none of the nearest neighbor is included, then check next-nearest neighbor to grow the module
      int size1 = NL.size();
      if(size1==size0)  {
	for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	  int j = (*p); 
	  
	  for(Nbl_itr q = A[j].begin(); q!= A[j].end(); q++) {
	    int k = (*q); 
	    if(!find(NL, k)) {
	      NL.push_back(k);
	      
	      double S0new, Snew;
	      CalModuleScore(NL, S0new, Snew); 
		
	      if(Snew>S*(1+r)) { 
		S = Snew;
		//if(S0new>S0*(1+r)) {
		//S0 = S0new;
		int j0 = GetHighestscoringCommonNeighbors(i, k);
		if(!find(NL,j0)) {
		  NL.push_back(j0);
		  CalModuleScore(NL, S0, S);
		}
	      }
	      else  
		NL.remove(k);
	    }
	  }
	}
      } // end of checking next-nearest neighbors of node i

      visited[i] = true;
    }// end of checking local enviorment of node i

  }// end of growing module with seed gene I


  // sort the node list, so that we can easily check if it has been obtained before.
  NL.sort();

}



/*
  Apparently, we can do this greedy search in parallel!
  But how!!!
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void State::HMS_Jia(double r, char* file)  
{
  cout << "\nSearch for the high-scoring modules using Jia's greedy algorithm: " << endl;

  vector<Nbl> AllNodelists;
  Modules.clear();

  double Smax = -1;
  for(int I=0; I<N; I++) {
    //cout << "Choose node " <<  I << " as the seed gene.\r";
    
      list<int> NL;
      NL.push_back(I); 
        
      ////////////////////////////////////////
      GreedyGrowth(NL, r);  
      ////////////////////////////////////////
      
      double s0, s;
      if(!find(AllNodelists, NL)) {
	AllNodelists.push_back(NL);
	
	CalModuleScore(NL, s0, s);
	
	module Module;
	Module.size= NL.size();
	Module.score0 = s0;
	Module.score = s;
	Module.size = NL.size();
	Module.nodelist = NL;
	
	Modules.push_back(Module);
	
	if(Module.score > Smax) {
	  Smax = Module.score;
	  cout << I << ' ' << Module.size << ' ' << Module.score << endl;
	}
      }
    
  }// end of growing modules for all nodes

  CalOverlapping();
  WriteGML(file);

  SaveModules(file);

  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Indeed, Jia's method yields much higher modules than Ideker's method. I think this is because in Jia's local greedy 
   algorithm, Jia makes full use of the network topology, while in Ideker's method, Ideker starts from a random initial state, 
   which requires zero biological knowledge. 

   However, Jia's method yields many overlapping modules. How shall we deal with this? (Note that Ideker's method will never generate 
   overlapping modules, but the score of those modules is much lower than Jia's method.) 

   We can do the following:
   (1) simply calculate the edge weight, i.e., counting the times that an edge appear in all the modules.   (2) apply frequent itemset mining technique to mine the frequent node sets (then check whether there are connected or not). Note that one can also perform the frequent sungraph mining (but there are so many different algorithms there, and I don't know which one is the best)
*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void State::CalOverlapping()
{

  int Nm = Modules.size(); 
  for (int a=0; a<Nm; a++) {
    vector<int> X; 
    for(Nbl_itr p = Modules[a].nodelist.begin(); p!= Modules[a].nodelist.end(); p++)
      X.push_back(*p);

    int n = X.size();
    for(int i=0; i<n; i++) {
      for(int j=i+1; j<n; j++) {
	if(find(A[X[i]],X[j])) {
	  LINK[edgeindex(X[i],X[j])].AddWeight(1.0); 
	  //cout << X[i] << ' ' << X[j] << ' ' <<  L[edgeindex(X[i],X[j])].GetWeight() << endl;
	}
      }
    }
    /*
      for(Nbl_itr p = Modules[a].nodelist.begin(); p!= Modules[a].nodelist.end(); p++) {
      int i = *p;
      for(Nbl_itr q = Modules[a].nodelist.begin(); q!= Modules[a].nodelist.end(); q++) {
      int j = *q;
      if(j!=i && find(A[i],j)) {
      L[edgeindex(i,j)].AddWeight(1.0); 
      cout << i << ' ' << j << ' ' <<  L[edgeindex(i,j)].GetWeight() << endl;
      }
      }
      }*/
  }
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
// check the connectivity of a node vector
bool State::CheckFNSConnectivity(vector<int>& NV)
{
  // We start from source node to do a BFS, if we can find all other nodes through nodes within NV, 
  // then this NV is connected. Otherwise it is not.  
  int source = NV[0]; 

  // assign each node a status variable: status[i] =1/0 : node i is in the vector NV or not  
  vector<char> status(N, 0);
  int n = NV.size();
  for(int a=0; a<n; a++) 
    status[NV[a]] = 1; 
  
  vector<bool> visited(N, false);

  queue<int> Q; // a FIFO queue   
  Q.push(source); //insert source into the queue
  while(!Q.empty())  {
    int u = Q.front(); 
    Q.pop();           
    visited[u] = true;
    //cout << NodeName[u] << ';'; // debug
    // loop over node u's adjacent nodes
    for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
      int v = (*p);
      // if this node is in the vector NV but has never been visited before.
      if(status[v]==1 && !visited[v])
	Q.push(v);            // put it into the queue
    }// end of loop over node u's adjacent nodes
  }// end of while loop

  bool connectivity = true;
  for(int a=0; a<n; a++) {
    if(!visited[NV[a]]) {
      connectivity = false;
      break;
    }
  }
  
  return connectivity;
}
////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////
// read the edgelist file, get the network
// then read the frequent node set file (which is the output file of the fpgrowth program)
// finally check for each frequent node set (FNS), whether it forms a connected subgraph or NOT
// This is equivalent to the frequent subgraph mining problem! 
void State::CheckFNSConnectivity(char* elistfile, char* fnsfile)
{
  cout << "Read the frequent node set file:\n";
  ifstream fin3(fnsfile, ios_base::in);
  if(!fin3) {cout << "Cannot open " << fnsfile << " for read."; exit(0);}
  char fname[256];
  sprintf(fname, "%s.connectivity", fnsfile);
  ofstream fout(fname, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  string S;
  while (getline(fin3, S)) {
    Splitter split(S, " ");
    int k = (int)split.size();
    int connectivity = 0;
    if(k==2) { // then this module contains only one node
      fout << S << ' ' << 1 << ' ' << 1 << endl;
    }
    else if(k>2) {
      vector<int> NV;
      for(int a=0; a<k-1; a++) {
	NV.push_back(NodeMAP[split[a].c_str()]);
      }
      connectivity = CheckFNSConnectivity(NV);

      // how to check whether the node set form a connected subgraph?
      // just do a BFS (breadth-first-search): 
      fout <<         connectivity << ' ' <<          k-1 << ' ' <<         S           << endl;
      cout << Red  << connectivity << ' ' << Green << k-1 << ' ' << Blue << S << Normal << endl;
    }
  }
  fin3.close();
  fout.close();

}
///////////////////////////////////////////////////////////////////////////////
