/**************************************************************************
 * adj_matrix.h
 *
 **************************************************************************/

#ifndef ADJ_MATRIX_H
#define ADJ_MATRIX_H

#include <vector>
#include <set>
#include <string>
#include "common.h"
#include "cost.h"
#include "occur.h"
#if defined (VFLIB)
	#include "argraph.h"
	#include "argedit.h"
	#include "attribute.h"
#endif

#define HALFKEY 16
#define aNULLKEY 0x00007fff
#define aMAKEKEY(x, y) ( ( (x) << HALFKEY ) |  (y)  )
#define aFIRSTKEY(x)   ( ( (x) >> HALFKEY ) )
#define aSECONDKEY(x)  ( ( (x)  << HALFKEY   >>  HALFKEY  ) )

typedef int KEYTYPE ;

class AdjMatrix 
{
private:
	vector<GLTYPE> matrix;
	vector<GNSIZE> ** SAdjList;
	short matrixSize;
	short eSize;      //each graph can only have 32k nodes
	short totalSup;   //we only support 32k graph transactions
    short bgf;

#ifdef SAVERESULTS
    //map<int, int>  graphs;
	//vector<occur*> occs;
	//void print(int i, string f1, string f2, vector<int> & g_reindex);
#if (SAVERESULTS == 2)
	GRAPHS1  TIDs;
#endif
#endif

#if defined VFLIB
	ARGraph<Point, void> * myg;
#endif

	void getMaxNode( vector<int> &v);
	void getMaxNode2( vector<int> &v);
	void last2bitil(int & i, int & l);
    int  lastbit();

	inline void assignLabel(int i, int j, GLTYPE label){

#ifdef SAFEMODE
		int s = size();
		if( i < 0 || i >= s || j < 0 || j >=s ){
			cout << "i is " << i << " j is " <<  j << endl;
			error("index outof boundary");
		}
#endif
		if( i < j) swap(i, j);
		matrix[i*(i+1)/2+j] = label;
	}

	inline void pushLabel(GLTYPE label){
		matrix.push_back(label);
	}

public:

	AdjMatrix(){ 
		//ocsize = 0; 
		SAdjList = NULL; 
		matrixSize=0;
		eSize = 0;
		totalSup = bgf = 0;
#ifdef VFLIB
  	    myg = NULL;
#endif
	}

	//constructor
	AdjMatrix(AdjMatrix * M, int index, GLTYPE el, GLTYPE nl){

		//int s = M->size();
		matrixSize = M->matrixSize;
		//eSize = M->eSize;
		//copy(M->matrix.begin(), M->matrix.end(), matrix.begin());
		matrix = M->matrix;
		if( nl != EMPTY_EDGE ) addNode(nl);
		addEdge(matrixSize-1, index, el);
		//ocsize =0;
		eSize = 0;
		totalSup = bgf = 0;

		SAdjList = NULL;
        buildAdjList(); 

#ifdef  VFLIB 
		myg = NULL;
#endif
	}
	
	AdjMatrix(AdjMatrix * M, int i, int j, GLTYPE el, GLTYPE nl){

		matrixSize = M->matrixSize;
		//eSize = M->eSize;
		matrix = M->matrix;
		if( i >=  matrixSize ) addNode(nl);
		addEdge(i, j, el);
		eSize = 0;
		totalSup = bgf = 0;
		SAdjList = NULL;
        buildAdjList(); 

#ifdef  VFLIB 
		myg = NULL;
#endif
	}

	//constructor
	AdjMatrix(AdjMatrix * M, Cost & lc){

		int s = M->size();
		matrixSize = M->matrixSize;
		//copy(M->matrix.begin(), M->matrix.end(), matrix.begin());
		matrix = M->matrix;
		for( s=0; s< lc.size(); s++)
			pushLabel(lc.node(s));

		totalSup = bgf = 0;
		SAdjList = NULL;
        buildAdjList(); 

#ifdef  VFLIB 
		myg = NULL;
#endif
	}

	~AdjMatrix() { 
	
		int i;
		if( SAdjList){ 
			for( i=0; i< size(); i++)   delete SAdjList[i];
            delete []SAdjList;
		}

#if (SAVERESULTS == 3)
		for(i=0; i< occs.size(); i++){
			delete occs[i];
		}
#endif

	}

	vector<GNSIZE> * getNeighbors(int index);
	void getNeighbors(int index, vector<GNSIZE> &v);

	int  isCanonicalForm1(); 
	int  isCanonicalForm2() ;
	int  isCanonicalForm();
    void setTotalSup (int s){ totalSup = s; }
	void setBackgroundF(int tc){ bgf = tc; }

    bool isPath();
	bool isTree(){ return ( matrixSize == (esize()+1) ); }
	inline GLTYPE getLabel(int i, int j) { 
#ifdef SAFEMODE
		if( i< 0 || i >= size() || j < 0 || j>= size()){ 
		  //print();
			cout << " i " << i << " j " << j << " size " << size() << endl;
			error("out of boundary");
		}
#endif
		if( i < j) swap(i,j);
		return matrix[i*(i+1)/2+j];
	}
	
	inline int  size() { return matrixSize; }
	int  esize();
	void addNode(char label);
	void addEdge(int i, int j, char label);
	void removeEdge(int i, int j);
	string canonicalForm();
	void buildAdjList();
	bool isSubgraphOf(AdjMatrix * m, int f);
	bool isIsomorphicOf(AdjMatrix *m){ return isSubgraphOf(m, 0) && m->isSubgraphOf(this,0) ; }
	void print(int i, string f1, string f2, string f3, map<int,int> & g_reindex);
    void print(int i, string f1, string f2);
	void print();
    bool denseEnough(int d, bool final);

#ifdef SAVERESULTS
    //void setGraphs(vector<occur*> & occs1);
#if (SAVERESULTS == 2)
	void setTIDs(GRAPHS1 & tids) { TIDs = tids; } 
#endif
	void setTIDs(GRAPHS1 & g,  int tsize);
#endif

#ifdef VFLIB 
	void setGraphForMatch();
	ARGraph<Point, void> * getARGraph(){ return myg; }
	bool isSubgraph(AdjMatrix * m);
	bool isIsomorphic(AdjMatrix * m);
	bool static my_visitor(int n, node_id ni1[], node_id ni2[], void *usr_data);
	int  findAllIsomophisms(AdjMatrix * m);
#endif

};

#endif
