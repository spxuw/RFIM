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

//#include "Histogram.h"

using namespace std;

typedef list<int> Nbl; // neighbor list
typedef Nbl::iterator Nbl_itr;
typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;
typedef vector<MATRIX> TENSOR;
const string numbers="0123456789";

//const string red="#ff0000";
//const string green="#00ff00";
//const string blue="#0000ff";




const int  INFINITE = INT_MAX;
const int  NIL      = -1;//INT_MAX;
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  vector<Link>  L;
  map<string,int> MAP; // the map between (i,j) and the edge index e

  int N;
  int E;
  vector<int> K;
  vector<Nbl> A;      // Adjacency list
  MATRIX AM;      // Adjacency Matrix

  int Ncopy;
  int Ecopy;
  vector<int> Kcopy;
  vector<Nbl> Acopy;  

  
  Rand rand;
  /////////////////////////////////////////////////////////////////////////////////////////////////
  list<int> LClist; // the vertex list of the largest component

  // To build a mapping between indices: we use the following two vertors
  // OriginalIndexof[i] give the original index of the largest component's i-th node 
  vector<int>  OriginalIndexofLCCnode;
 
  // LCCIndexof[i] give the in-largest-component index of the original network's i-th node
  vector<int>  LCCIndexofOriginalnode; // if IndexinLC[i] = -1, this means it doesn't belong to the lc
  // otherwise, its value just gives the index

  int Ncc; // # of connected components
  vector<Nbl>  AllComponents; // connected component
  vector<int>  ccindex; // which connected component each node belongs to 

  vector<bool> LC; // whether the node belongs to the largest component
  int          lcc_index; // the index of the largest component
  int          Nlc; // number of nodes in the largest connected component
  double       Nac; // average size of connected components

  double       nlc; // :=Nlc/N
  int          Elc; // number of edges in the largest connected component
  double       mlc; // edge density for the largest connected component 
  /////////////////////////////////////////////////////////////////////////////////////////////////


  MATRIX M; // neighborhood matrix
  int diameter;

};

inline int Graph::edgeindex(int i, int j)
{
  stringstream sst; sst << i << ">" << j;
  return MAP[sst.str()];
}



#endif
