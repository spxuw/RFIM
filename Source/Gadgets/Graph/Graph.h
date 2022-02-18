#if !defined (GRAPH_H)
#define GRAHP_H


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <queue> 
#include <deque>
#include <map>
#include <stack>
#include <set>

#include <limits.h>
#include <stack> 
#include "Link.h"
#include "find.h"
#include "Rand.h"
#include "color.h"



using namespace std;

typedef list<int> Nbl; 
typedef Nbl::iterator Nbl_itr;
typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;
typedef vector<MATRIX> TENSOR;
const string numbers="0123456789";








const int  INFINITE = INT_MAX;
const int  NIL      = -1;
const int  INF      = INT_MAX;


class Graph
{
 public:
  Graph();
  Graph(MATRIX& LinkWeight, double t);
  ~Graph();

  void Erdos_Renyi(int n, double p, int seed);

  void Backup();
  void Restore();

  void AddLink(int i, int j);
  int  edgeindex(int i, int j);

  void RemoveLink(int i, int j);
  void AddBackLink(int i, int j);
  void RemoveNode(int i);
  void Remove(set<int>& V);

  bool Plot(char* fname);  

  void DFS();
  void DFS(int u, vector<bool>& visited, Nbl&  Component);

  void BFS(int source, VECTOR& d);
  int  BFS(int source, int target);
  void CalNeighborhoodMatrix();
  double GetNetworkDistance(Graph& G0);
  double GetMatrixDistance(Graph& G0);

  
  vector<Link>  L;
  map<string,int> MAP; 

  int N;
  int E;
  vector<int> K;
  vector<Nbl> A;      
  MATRIX AM;      

  int Ncopy;
  int Ecopy;
  vector<int> Kcopy;
  vector<Nbl> Acopy;  

  
  Rand rand;
  
  list<int> LClist; 

  
  
  vector<int>  OriginalIndexofLCCnode;
 
  
  vector<int>  LCCIndexofOriginalnode; 
  

  int Ncc; 
  vector<Nbl>  AllComponents; 
  vector<int>  ccindex; 

  vector<bool> LC; 
  int          lcc_index; 
  int          Nlc; 
  double       Nac; 

  double       nlc; 
  int          Elc; 
  double       mlc; 
  


  MATRIX M; 
  int diameter;

};

inline int Graph::edgeindex(int i, int j)
{
  stringstream sst; sst << i << ">" << j;
  return MAP[sst.str()];
}



#endif
