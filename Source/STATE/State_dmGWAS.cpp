#include "State.h"

using namespace std;






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






void State::GreedyGrowth(Nbl& NL, double r)
{
  double S0, S, S0new, Snew;
  CalModuleScore(NL, S0, S);

  vector<bool> visited(N, false);
  for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) {
    int i = (*p);
    if(!visited[i]) {
      
      int size0 = NL.size();
      
      for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
	int j = (*p); 
	
	if(!find(NL, j)) {
	  NL.push_back(j);
	  
	  CalModuleScore(NL, S0new, Snew); 
	  
	  if(Snew>S*(1+r))
	    S = Snew;
	  
	  
	  else  
	    NL.remove(j); 
	}
      }
      
      
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
      } 

      visited[i] = true;
    }

  }


  
  NL.sort();

}





void State::HMS_Jia(double r, char* file)  
{
  cout << "\nSearch for the high-scoring modules using Jia's greedy algorithm: " << endl;

  vector<Nbl> AllNodelists;
  Modules.clear();

  double Smax = -1;
  for(int I=0; I<N; I++) {
    
    
      list<int> NL;
      NL.push_back(I); 
        
      
      GreedyGrowth(NL, r);  
      
      
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
    
  }

  CalOverlapping();
  WriteGML(file);

  SaveModules(file);

  
}



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
	  
	}
      }
    }
   
  }
 
}





bool State::CheckFNSConnectivity(vector<int>& NV)
{
  
  
  int source = NV[0]; 

  
  vector<char> status(N, 0);
  int n = NV.size();
  for(int a=0; a<n; a++) 
    status[NV[a]] = 1; 
  
  vector<bool> visited(N, false);

  queue<int> Q; 
  Q.push(source); 
  while(!Q.empty())  {
    int u = Q.front(); 
    Q.pop();           
    visited[u] = true;
    
    
    for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
      int v = (*p);
      
      if(status[v]==1 && !visited[v])
	Q.push(v);            
    }
  }

  bool connectivity = true;
  for(int a=0; a<n; a++) {
    if(!visited[NV[a]]) {
      connectivity = false;
      break;
    }
  }
  
  return connectivity;
}











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
    if(k==2) { 
      fout << S << ' ' << 1 << ' ' << 1 << endl;
    }
    else if(k>2) {
      vector<int> NV;
      for(int a=0; a<k-1; a++) {
	NV.push_back(NodeMAP[split[a].c_str()]);
      }
      connectivity = CheckFNSConnectivity(NV);

      
      
      fout <<         connectivity << ' ' <<          k-1 << ' ' <<         S           << endl;
      cout << Red  << connectivity << ' ' << Green << k-1 << ' ' << Blue << S << Normal << endl;
    }
  }
  fin3.close();
  fout.close();

}

