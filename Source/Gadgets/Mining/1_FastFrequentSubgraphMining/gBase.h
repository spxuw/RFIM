#ifndef GBASE_H
#define GBASE_H

#include <vector>
#include <map>
#include "adj_matrix.h"
#include "occur.h"

class gBase
{
private:
	
	//inner class for input processing
	class Edge {
	public:
		int  node1_no; // label of node1 < label of node2
		int  node2_no;
		char label;
		int  key;

		Edge(int n1, int n2, char e, char L1, char L2) : node1_no(n1), node2_no(n2) 
		{ label = e; key = L2 << 16 | label << 8 | L1; }   // L1 < L2
	};


	//date strutures for graph database and input processing
	vector<AdjMatrix*> gb;
        string nodefile;
        string edgefile;
	int threshold;
	map<int, set<int> *> edges_freq;
	map<int, vector<occur *> *, greater<int> > edges_occur;
	map<int, vector<occur *> *, greater<int> > node_occur;
	vector< vector<int> *> v_reindex;
	vector< int > g_reindex;

	void scanNode(ifstream &fnode, vector< vector<char> *> &ng, map<char, set<int>*> &nf);
	void scanEdge(ifstream &fedge, vector< vector<char> *> &ng, vector< vector<Edge *>* > &eg);
	void scanII(vector< vector<char> *> &n, vector< vector<Edge *>* > &e, map<char, set<int>*> &nf);
	void scanIII();
	void swap(int *x, int *y);
	void printEdge(Edge *e);
	void printKey(int k);
	
public:
	gBase(string &f1, string &f2, int f): nodefile(f1), edgefile(f2), threshold(f) {scan(); }
	~gBase();

	void setAdditionalInfo(string &f);
        void scan();
	char edge (int gi, int n1, int n2) {return gb[gi]->getLabel(n1, n2); }
	AdjMatrix *graph (int gi) { return gb[gi]; }
	map<int, vector<occur *>*, greater<int> > &getMap() { return edges_occur; }
	map<int, vector<occur *>*, greater<int> > &getNodeMap() { return node_occur; }
	
	void deleteMap(){ 
		map<int, vector<occur *> *, greater<int> >::iterator ip = edges_occur.begin();
		for(; ip != edges_occur.end(); ip++){
			vector<occur*> * vc = ip->second;
			for( int i=0; i< vc->size(); i++) delete (*vc)[i];
			delete vc;
		}
	}
    
	void deleteNodeMap(){
		map<int, vector<occur *> *, greater<int> >::iterator ip = node_occur.begin();
		for(; ip != node_occur.end(); ip++){
			vector<occur*> * vc = ip->second;
			for( int i=0; i< vc->size(); i++) delete (*vc)[i];
			delete vc;
		} 
	}
	void print();

	//add by luke
	void getNeighbors( int gi, int n1, vector<GNSIZE> &v) { gb[gi]->getNeighbors(n1, v); }
	//void getNeighbors( int gi, int n1, vector<GNSIZe> &v, set<int> & filter) { gb[gi]->getNeighbors(n1, v, filter); }
	int  getThreshold(){  return threshold; }
	void removeEdge(int key);
	int  size() { return gb.size(); }
	vector<int> & getGIndex(){ return g_reindex; }
	void getGIndex(map<int, int>& gMap){ 
		for(int i =0; i< g_reindex.size(); i++)
			gMap[g_reindex[i]] = i; 
	}
	bool isSubgraph(int i, AdjMatrix * m ) { return gb[i]->isSubgraphOf(m, 0); } 
};

class Comp
{
public:
        bool operator()(  occur * const & oc1,  occur *  const & oc2);
};

#endif

