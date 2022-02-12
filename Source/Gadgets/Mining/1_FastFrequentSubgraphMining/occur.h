#ifndef OCCUR
#define OCCUR

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "common.h"

//#define  SETSTAGE(x)   ((x) | 1)
//#define  CLEARSTAGE(x) ((x) & 0xffffff00)
//#define  GETSTAGE(x)   ((x) & 0x000000ff) 

class occur
{
private:

	vector<GNSIZE>		nodes;
	int					graph;
	
public:

	//constrcutor
	occur( vector<GNSIZE> & n1, int g1){
		nodes = n1;
		graph = g1;	
	}

	//deconstrcutor
	~occur(){
		;
	}

	inline int getGraph(void) { return graph ; }
	inline int getNodes(int i){ return nodes[i]; }
	inline int size() { return nodes.size(); }
	
	void print(){
		int i;
		for(i = 0; i<nodes.size(); i++){
			cout << nodes[i] << "   "; 
			if( (i+1) % 5 == 0) cout<< endl;
		}
		if( i % 5 != 0 ) cout << endl;

		cout << " graph is " << graph  << endl;
	}
};

#endif

