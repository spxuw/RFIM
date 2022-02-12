#ifndef ULLMAN_H
#define ULLMAN_H

#include <vector>
#include <set>
#include <string>
#include "occur.h"
#include "common.h"
#include "adj_matrix.h"
#include "pattChild.h"

#ifdef INDUCEDS
	typedef map<string, pattern*> CHILDREN;
#else
	typedef map<int, pattern*>	  CHILDREN;
#endif

typedef vector<occur*> OCCS ;
class Umapper{
private:
	BYTE *M;
	int n1, n2;

public:

	//constructor
	Umapper(int s1, int s2){
		n1 = s1; n2 = s2;
		M = new BYTE[n1*n2];
		memset(M, UNEDGE, n1*n2);
	}
	
	Umapper(Umapper * M1){
		n1 = M1->n1; n2 = M1->n2; 
		M = new BYTE[n1*n2];
		memcpy(M, M1->M, n1*n2);
	}

	~Umapper(){ delete [] M; }

	inline void setEdge(int i, int j){ M[i*n2+j] = UEDGE; }
	inline void clearEdge(int i, int j){ M[i*n2+j] = UNEDGE; }
	void clearEdges(int i, int j){
		for(int k=0; k< n2; k++) {
			if( k!= j)	clearEdge(i, k);
		}
		//M[i*n2+k] = UNEDGE; 
	}
	
	inline BYTE getEdge(int i, int j){ return M[i*n2+j];}
	inline bool hasEdge(int i, int j){ return getEdge(i, j) == UEDGE ; }

};

class UllmanIsomorphism
{
private:

	AdjMatrix *g1, *g2;
	pattern * pat;
	int n1, n2;
	int gi;

	void Ullman_initial(Umapper * M0, vector<int> &Col, vector<int> & F);
	bool Ullman_backtrack(Umapper * M0, int d, vector<int> &Col, vector<int> & F);
//	void proposeChild(vector<int> &Col, vector<int> & F);

/*
#ifdef INDUCEDS
	void proposeIndu(vector<int> &Col, vector<int> & F);
#else
	void proposeInner(vector<int> &Col, vector<int> & F);
	void proposeOuter(vector<int> &Col, vector<int> & F);
#endif
*/

public:
	
	//constructor
	UllmanIsomorphism(pattern * p1, AdjMatrix * s2, int gi){
		pat = p1; //g1 = pat->getAdjMatrix();
		g2 = s2; this->gi = gi;
		n1 = g1->size(); n2 = g2->size();
	}

	UllmanIsomorphism(AdjMatrix * s1, AdjMatrix * s2){
		pat = NULL; g1 = s1;
		g2 = s2; this->gi = 0;
		n1 = g1->size(); n2 = g2->size();
	}

	//deconstructore
	~UllmanIsomorphism(){};

	bool Ullman_isomorphism( );


};

#endif

