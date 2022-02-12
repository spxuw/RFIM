#include "State.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////
void State::WriteGML(char* fname)
{
  char filename[256];
  sprintf(filename, "%s.gml", fname);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  //time(&rawtime);

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; //so that any number will be printed in this fixed notation, because Cytoscape cannot read scientific notation, e.g., 1e-3
  fout << setprecision(15); //otherwise, 1e-7 would be 0.000000

 
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
    int i = LINK[e].GetHead(); // starting node of the l-th link 
    int j = LINK[e].GetTail(); // ending node of the l-th link 
    fout << "edge [" << endl;
    fout << "source " << i << endl;
    fout << "target " << j << endl;
    fout << "EdgeWeight " <<  LINK[e].GetWeight() << endl;
    fout << "] " << endl;
  }
  fout << "]";
  fout.close();

}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// loop over node u's adjacent nodes to check if there is any neighbor=DIR
bool State::CheckNeighbors(char dir, int u) 
{
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    if(spin[*p]==dir)
      return true;
  }
  
  return false; 
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// write only those components with spin values= dir//, and benchmark genes highlighted
void State::WriteGML(char dir, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d.gml", fname, (int)dir);
  ofstream fout(filename, ios_base::out);

  sprintf(filename, "%s-spin%d.txt", fname, (int)dir);
  ofstream fout2(filename, ios_base::out);


  time_t rawtime;
  //time(&rawtime);

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; //so that any number will be printed in this fixed notation, because Cytoscape cannot read scientific notation, e.g., 1e-3

  vector<bool> highlighted(N,false); // those nodes will be highlighted in the PPI

  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); // starting node of the l-th link 
    int j = LINK[e].GetTail(); // ending node of the l-th link 
    //   if(spin[i]==dir && spin[j]==dir) {
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
    //if(spin[i]==dir || CheckNeighbors(dir, i)) {
    if(highlighted[i]) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      //fout << "benchmark "<< (int)benchmark[i] << endl; 
      fout << "] " << endl;
    }

    if(spin[i]==dir)
      fout2 << NodeName[i] << endl; 
  }
  

  fout << "]";
  fout.close();

  fout2.close();

}
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
// write only those components with spin values= dir and size > smin//, and benchmark genes highlighted
void State::WriteGML(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  //time(&rawtime);

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; //so that any number will be printed in this fixed notation, because Cytoscape cannot read scientific notation, e.g., 1e-3

  for(int i=0; i<N; i++) {
    if(spin[i]==dir && Modules[moduleindex[i]].size >smin) {
      fout << "node [" << endl;
      fout << "id " << i << endl;
      fout << "label \"" << NodeName[i] << "\"" << endl;
      fout << "K "<< K[i] << endl;    
      fout << "Spin "<< (int)spin[i] << endl; 
      fout << "NodeWeight "<< h[i] << endl; 
      fout << "NodePvalue "<< PVALUE[i] << endl; 
      //fout << "benchmark "<< (int)benchmark[i] << endl; 
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); // starting node of the l-th link 
    int j = LINK[e].GetTail(); // ending node of the l-th link 
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
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// write only those components with spin values= dir and size > smin//, and benchmark genes highlighted
void State::WriteGML_HightlightBenchmarkGenes(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d-highlight.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  //time(&rawtime);

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; //so that any number will be printed in this fixed notation, because Cytoscape cannot read scientific notation, e.g., 1e-3

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
      //fout << "benchmark "<< (int)benchmark[i] << endl; 
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); // starting node of the l-th link 
    int j = LINK[e].GetTail(); // ending node of the l-th link 
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
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// write only those components with spin values= dir and size > smin//, and benchmark genes highlighted
// also draw their first nearest neighbors (regardless of their spin values) YYL 10/15/2015
void State::WriteGML_HightlightBenchmarkGenes_NearestNeighbors(char dir, int smin, char* fname)
{
  char filename[256];
  sprintf(filename, "%s-spin%d-smin%d-highlight.gml", fname, (int)dir, smin);
  ofstream fout(filename, ios_base::out);

  time_t rawtime;
  //time(&rawtime);

  fout << "Creator \"YYL on " << ctime(&rawtime) << "\"" << endl;
  fout << "graph [" << endl;
  fout << "comment \"This is a graph with the edge weights and node weights highlighted\"" << endl;
  fout << "directed 0" << endl;

  fout << fixed; //so that any number will be printed in this fixed notation, because Cytoscape cannot read scientific notation, e.g., 1e-3


  // check if a node has a flipped neighbor
  vector<bool> HasFlippedNeighbors(N, false); 
  for(int i=0; i<N; i++) {
    for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) {
      int j = *p;
      //cout << i << ',' << j << endl; //test
      if((spin[j])==dir) {
	HasFlippedNeighbors[i] = true; 
      }
    }
    //cout << i << ' ' << (int)(HasFlippedNeighbors[i]) << endl; //test
  }
  //cout << "test\n ";

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
      //fout << "benchmark "<< (int)benchmark[i] << endl; 
      fout << "] " << endl;
    }
  }
  
  int Ne = LINK.size();
  for(int e=0; e<Ne; e++) {
    int i = LINK[e].GetHead(); // starting node of the l-th link 
    int j = LINK[e].GetTail(); // ending node of the l-th link 
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
/////////////////////////////////////////////////////////////////////////////////////////////////
