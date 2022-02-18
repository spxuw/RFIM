#include "State.h"

using namespace std;


void State::WriteGML(char* fname)
{
  char filename[256];
  sprintf(filename, "%s.gml", fname);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; 
  fout << setprecision(15); 

 
  for(int i=0; i<N; i++) {
   
    fout << "node [" << endl;
    fout << "id " << i << endl;
    fout << "label \"" << NodeName[i] << "\"" << endl;
    fout << "K "<< K[i] << endl;    
    fout << "Spin "<< (int)spin[i] << endl; 
    fout << "NodeWeight "<< h[i] << endl; 
    fout << "NodePvalue "<< PVALUE[i] << endl; 
    fout << "] " << endl;
  }
 
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); 
    int j = LINK[e].GetTail(); 
    fout << "edge [" << endl;
    fout << "source " << i << endl;
    fout << "target " << j << endl;
    fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
    fout << "] " << endl;
  }
  fout << "]";
  fout.close();

}







bool State::CheckNeighbors(char dir, int u) 
{
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    if(spin[*p]==dir)
      return true;
  }
  
  return false; 
}






void State::WriteGML(char dir, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d.gml", fname, (int)dir);
  ofstream fout(filename, ios_base::out);

  sprintf(filename, "%s-spin%d.txt", fname, (int)dir);
  ofstream fout2(filename, ios_base::out);


  time_t rawtime;
  

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; 

  vector<bool> highlighted(N,false); 

  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); 
    int j = LINK[e].GetTail(); 
    
    if(spin[i]==dir || spin[j]==dir) {
      fout << "edge [" << endl;
      fout << "source " << i << endl;
      fout << "target " << j << endl;
      fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
      fout << "] " << endl;

      highlighted[i]=true;
      highlighted[j]=true;
    }
  }


  for(int i=0; i<N; i++) {
    
    if(highlighted[i]) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      
      fout << "] " << endl;
    }

    if(spin[i]==dir)
      fout2 << NodeName[i] << endl; 
  }
  

  fout << "]";
  fout.close();

  fout2.close();

}








void State::WriteGML(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; 

  for(int i=0; i<N; i++) {
    if(spin[i]==dir && Modules[moduleindex[i]].size >smin) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); 
    int j = LINK[e].GetTail(); 
    if(spin[i]==dir && spin[j]==dir && Modules[moduleindex[i]].size >smin) {
      fout << "edge [" << endl;
      fout << "source " << i << endl;
      fout << "target " << j << endl;
      fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
      fout << "] " << endl;
    }
  }

  fout << "]";
  fout.close();

}





void State::WriteGML_HightlightBenchmarkGenes(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d-highlight.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; 

  for(int i=0; i<N; i++) {
    if(spin[i]==dir && Modules[moduleindex[i]].size >smin) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      fout << "BenchmarkGene "<< Color[i] << endl; 
      
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); 
    int j = LINK[e].GetTail(); 
    if(spin[i]==dir && spin[j]==dir && Modules[moduleindex[i]].size >smin) {
      fout << "edge [" << endl;
      fout << "source " << i << endl;
      fout << "target " << j << endl;
      fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
      fout << "] " << endl;
    }
  }

  fout << "]";
  fout.close();

  sprintf(filename, "%s-spin%d-smin%d-highlight.txt", fname, (int)dir, smin);
  ofstream fout2(filename, ios_base::out);

 for(int i=0; i<N; i++) {
    if(spin[i]==dir && Modules[moduleindex[i]].size >smin) {
      fout2 << NodeName[i] << ' ' << K[i] << ' ' << (int)spin[i] << ' ' << h[i] << ' ' << PVALUE[i] << ' ' << Color[i] << endl; 
    }
 }
 fout2.close();


}






void State::WriteGML_HightlightBenchmarkGenes_NearestNeighbors(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d-highlight.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; 


  
  vector<bool> HasFlippedNeighbors(N, false); 
  for(int i=0; i<N; i++) {
    for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
      int j = *p;
      
      if((spin[j])==dir) {
	HasFlippedNeighbors[i] = true; 
      }
    }
    
  }
  

  for(int i=0; i<N; i++) {
    if((spin[i]==dir && Modules[moduleindex[i]].size >smin) || HasFlippedNeighbors[i]) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      fout << "BenchmarkGene "<< Color[i] << endl; 
      
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); 
    int j = LINK[e].GetTail(); 
    if((spin[i]==dir && spin[j]==dir && Modules[moduleindex[i]].size >smin) || 
       (HasFlippedNeighbors[i] && HasFlippedNeighbors[j])
       ) {
      fout << "edge [" << endl;
      fout << "source " << i << endl;
      fout << "target " << j << endl;
      fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
      fout << "] " << endl;
    }
  }

  fout << "]";
  fout.close();

  sprintf(filename, "%s-spin%d-smin%d-highlight.txt", fname, (int)dir, smin);
  ofstream fout2(filename, ios_base::out);

 for(int i=0; i<N; i++) {
    if(spin[i]==dir && Modules[moduleindex[i]].size >smin) {
      fout2 << NodeName[i] << ' ' << K[i] << ' ' << (int)spin[i] << ' ' << h[i] << ' ' << PVALUE[i] << ' ' << Color[i] << endl; 
    }
 }
 fout2.close();


}

