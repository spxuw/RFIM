#pragma warning (disable:4786) 
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

using namespace std;   

#include "gBase.h"
#include "occur.h"
#include "common.h"
#include "Ullman.h"
 

void UllmanIsomorphism :: Ullman_initial(Umapper * M0, vector<int> &Col, vector<int> & F)
{
	//set up the M0
	int i;
	for(i=0; i< n1; i++){
		int t = g1->getNeighbors(i)->size();
		int nl1 = g1->getLabel(i, i);
		for(int j=0; j< n2; j++){
			if( t <= g2->getNeighbors(j)->size() &&
				nl1 == g2->getLabel(j, j)	) 
				M0->setEdge(i, j);
		}
	}
	//set Col and F
	for( i=0; i< n1; i++) F.push_back(-1);
	for( i=0; i< n2; i++) Col.push_back(0);
}

/*void UllmanIsomorphism:: refine (Umapper * M1, vector<int> & F, int d)
{
	for(int i=
}*/

bool UllmanIsomorphism:: Ullman_backtrack(Umapper * M0, int d,  
										  vector<int> &Col, vector<int> & F)
{
#ifdef DEBUG6
	cout << "d is n1 n2" << d << "  " << n1 << "  " << n2;
	for(int i=0; i< n1; i++){
		cout << F[i] << "  " ;
	}
	cout << endl;
#endif
	bool flag; 
	if( d >= n1){
		//proposeChild(Col, F);
		return true;
	}
	for(int j=0; j< n2; j++){
		if( M0->hasEdge(d, j) && Col[j]==0 &&
			g1->getLabel(d,d) == g2->getLabel(j, j) ){
			flag = true;
			int k;
			for( k=0; k< d; k++){
				//for inducded subgraph, empty edge also meaningful
#ifdef INDUCEDS
				if( g1->getLabel(d,k) != g2->getLabel(j, F[k]) ){
#else
				if( g1->getLabel(d, k) != EMPTY_EDGE && g1->getLabel(d,k) != g2->getLabel(j, F[k]) ){
#endif
					flag = false; break;
				}  
			}

			if( flag ) { //match so far 
				//vector<int> Col1 = Col;	
				//vector<int> F1   = F;		
				Col[j] = 1; F[d] =j;
				Umapper * M1 = new Umapper(M0);
				M1->clearEdges(d, j);

				//refine the mapper 
				for(int ii=d+1; ii< n1; ii++){
					if( g1->getLabel(ii, d) != EMPTY_EDGE)
						for(int jj=0; jj< n2; jj++){
							if( M1->hasEdge(ii, jj) && g1->getLabel(ii, d) != g2->getLabel(jj, j) )
								M1->clearEdge(ii, jj);
						}  
				}

				bool result = Ullman_backtrack(M1, d+1, Col, F);
				Col[j] = 0; F[d] =-1;
				delete M1;
				if( result ) return true;
			}
		}
	}
	return false;
}

bool UllmanIsomorphism :: Ullman_isomorphism()
{
//	cout << "start isomorphism search" << endl;

	/*
	cout << "m1" << endl;
	g1->print();

	cout << "m2 " << endl;
	g2->print();
	*/

        Umapper * M0 = new Umapper(n1, n2);
	vector<int> Col;
	vector<int> F;
	Ullman_initial(M0, Col, F);
  //      cout << "start backtrack" << endl;
	bool f = Ullman_backtrack(M0, 0, Col, F);
	delete M0;

	//cout << "result is " << f << endl;
	return f;
}

//#endif
