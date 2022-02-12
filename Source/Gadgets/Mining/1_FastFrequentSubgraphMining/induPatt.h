#ifndef INDUPATT_H
#define INDUPATT_H

#include <set>
#include <map>
#include <vector>

#include "common.h"
#include "adj_matrix.h"
#include "cost1.h"
#include "cost.h"
#include "occur.h"
#include "gBase.h"

//#define iMAKEKEY(x, y) ( ( ( (DLONG)(x) ) << HALFDLONG ) |  ((DLONG)(y) ) )
//#define iSECONDKEY(x)  ( (int) ( (x)  << HALFDLONG  >>  HALFDLONG ) )

typedef map<int, vector<occur*> *, greater<int> > EDGES;
typedef vector<occur *>  OCCS;
//typedef set<int> GRAPHS;

class iPattern
{
private:
	AdjMatrix			* M;
	COCCS				* coccs;
	GRAPHS1				* gids;
	iPattern			* par;
	gBase				* gb;
	MyTimer				& Timers;
	Register			& Registers;
	int					optimal;
	short				stage;
	short				size;
	Cost				mycost;
	vector<iPattern *>	child;

	//used for initialization
	void setChild(EDGES & cm, EDGES & nm);
	iPattern * setUpChild(vector<int> u, EDGES & cm, EDGES & gnodes, int nlabel);
	void setUpEdges(iPattern * p2, OCCS * es, TRA & tra);
	void basicSetup(iPattern * p1) { 
		M = p1->M;	coccs = new COCCS();	gids = new GRAPHS1(); 
		par =p1;	stage = 0;	gb = p1->gb; optimal = 0; 
	}

	//used by join two patterns
	void join(iPattern * p2, int freq);
	void propose_outer2(iPattern * p2, int freq);
	void proposeOuterChild(int f);
    void scan(map<int, iPattern*, greater<int> > & cand, int f, vector<iPattern*> & result);

	//functinos for find all isomorphism
#ifdef GSPANFLAG
	void findAllIsomorphisms();
#endif

	//others 
	void constructM (int k );
	void fixM(void);
	int  isCarnonicalForm(void);
	
	//some interface between patterns
	inline int  sup() const { return gids->size(); } 
	inline int  getSize() { return size; }
	void sprint(int level = 0);
	void fsm2( vector<AdjMatrix * > & result, int freq, int level, int & total);

public:

	//constructor
	iPattern(gBase *gb1, MyTimer & t, Register & mr, EDGES & cm, EDGES & nm) : gb(gb1), Timers(t) , 
		Registers(mr) 
	{
		M = NULL;				coccs = new COCCS();			gids = new GRAPHS1(); 
		par = NULL;				stage = 0;						size = 0;				
		optimal = 0;			setChild(cm, nm);			
	}
	
	iPattern( iPattern * p1, int k, int s1):Timers(p1->Timers), Registers(p1->Registers)  {
		size = s1;		basicSetup(p1);	
		if( size == 1)  mycost.append(COSTNL(k));
		else{
			mycost.append(COSTEL(k));		mycost.append(COSTNL(k));	
		}
		constructM(k);
	}

	iPattern( iPattern * p1,  Cost & lc, int s1 ):Timers(p1->Timers), Registers(p1->Registers){
		size = s1;	basicSetup(p1);		mycost = lc;
	}

	//deconstructor
	~iPattern(){
		if( stage == 1 || stage == 3 || stage == 4){
	       		delete M;
		}
	
		if( stage > 4) error("stage is too large");
		delete gids;	
		delete coccs;
	}

	//member functions
	void fsm1( vector<AdjMatrix * > & result, int freq, int level);
	void print();

};

typedef vector<iPattern *> iCHLTYPE ;
typedef	map<int, iPattern * , greater<int> > CANDIDATE;

#endif
