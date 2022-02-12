#ifndef PATCHILD_H
#define PATCHILD_H

#include "common.h"
#include "adj_matrix.h"
#include "cost1.h"
#include "occur.h"
#include "gBase.h"

#define TURNON(x, y) ((x)[(y) /32]) |=  ( 1 << (31- ( (y) % 32) ) )
#define GETTID(x, y) ((x)[(y) /32]) &=  ( 1 << (31- ( (y) % 32) ) )

typedef map<int, vector<occur*> *, greater<int> > EDGES;
typedef vector<occur *>  OCCS;

#define SCANOCCURRENCES 	int edgeSize = M->esize(), j;  \
	pattern * currp = this; \
	vector<int> gindices(edgeSize, INULLINDEX); \
	vector<short> checkpoint(edgeSize+1, 0); \
	vector<COCCS*> coccslist; \
	int fkey, skey, fkey1, gi, level, i; \
	for( i=0; i<= edgeSize; i++){ \
	if( THIRDKEY((*(currp->coccs))[0]) != GNULLINDEX ){ \
	    checkpoint[i] = 1; }\
		coccslist.push_back(currp->coccs); \
		currp = currp->par; \
	} \
	vector<GNSIZE> filter(size, 0); \
	COCCS :: iterator ip = coccs->begin(); \
	for( i=0; ip != coccs->end(); i++ , ip++){\
		fkey = FIRSTKEY((*coccs)[i]), skey = THIRDKEY((*coccs)[i]); \
		gi = SECONDKEY((*coccs)[i]); \
		level = size-1; if( checkpoint[0] ){	filter[level--]= skey; } \
		for( j=0; j< edgeSize && fkey != gindices[j]; j++){ \
			gindices[j] = fkey;  \
			fkey1 = FIRSTKEY((*(coccslist[j+1]))[fkey]);  \
			if( checkpoint[j+1]){   \
				filter[level--]= THIRDKEY((*(coccslist[j+1]))[fkey]); } \
			fkey = fkey1; \
		} 

class pattern
{
private:
	AdjMatrix			* M;
	pattern				* par;
	char				optimal;
	char				stage;
	short				size;
	cost1				mycost;
	static gBase		* gb;
	static MyTimer		* Timers;
	static Register		* Registers;
	vector<pattern *>	child;
	vector<pattern *>	newchild;
	GRAPHS1				gids;
	COCCS				* coccs;
	static DLONG		totalPatterns;
		
#ifdef MAXIMALPATTERN 
	list<int>          edgeCan;
	list<GLTYPE>       edgeLabels;
#endif

	//funcions for proposing child
	void fsm2( vector<AdjMatrix * > & result, int freq, int level, 
		int & total, int density, string nodef, string edgef, 
		string header, int sizelimit, int sizeup, map<int,int> &gIndex);
	bool proposeOuterChild(int f, vector<pattern*> & outChild);
    void scan(map<int, pattern*, greater<int> > & cand, int f, vector<pattern*> & result);
	void fixTRA();

	//used for initialization
	pattern * setUpChild(vector<int> u, EDGES & cm, EDGES & gnodes, int nlabel);
	void setUpEdges(pattern * p2, OCCS * es, TRA & tra);
	
	//used by join two patterns
	void join(pattern * p2, int freq);
	bool join1(pattern * p2, pattern * p3, int type, int freq, bool mf);
	bool propose_inner(pattern * p2, pattern * p3, int freq);
	bool propose_inner_outer(pattern * p2, pattern * p3, int freq);
	bool propose_outer2(pattern * p2, pattern * p3, int freq);
	void populate_child( pattern * p3, vector<int> * s3, int t, occur * oc1);

#ifdef MAXIMALPATTERN
	void fillMoreEdge(vector<GNSIZE> & nodes, AdjMatrix * gbi);  //used by Maximal pattern mining
	bool isInnerMaximal();
	void initialEdgeCan();
	bool isMaximal(int threshold);
	bool countSuperGraph(vector<GNSIZE> & nodes, int gi, CGOCC & counter, int threshold);
#endif

#if    (TREEBASEDSCAN == 1)
	bool isTreeClosed();
	bool countSuperTree(vector<GNSIZE> & nodes);
#endif

	//functinos for find all isomorphism
#ifdef GSPANFLAG
	void findAllIsomorphisms();
#endif

	//save results
#ifdef SAVERESULTS
	void recordOccs(vector<occur*> & loccs);
#endif

	//others 
	void constructM (void);
	void fixM(void);
	int  isCarnonicalForm(void);
    bool denseEnough(int density, int sizeup, bool f);	
	//some interface between patterns
	inline int  sup() const { /*return graphs->size();*/ return gids.size(); } 
	inline int  getSize() { return size; }
	inline int  getCost() { return mycost.intValue(); }
	void sprint(int level = 0);

	//related Tree based scanning
	void frequentElement(vector<pattern*> & result, int freq);
	void enuFreqEle(FREELE & candidates, vector<pattern*> & results, int freq);
	void scanElements(FREELE & candidates, int freq);
	void addInstances(ELEFRENQC & counter, vector<GNSIZE> & occ, int gi, int instanceID );
	int  getSupport( IIDTYPE * iid);

public:

	//constructor
	pattern()
	{
		M = NULL;			coccs = new COCCS();			par = NULL;			
		stage = 0;			size = 0;			   		    optimal = 0;				
	}
	
	pattern( pattern * p1, int k, int s1) {
		M = p1->M;			coccs = new COCCS();			
		par =p1;			mycost.setCost(k);				stage = 0;						
		size = s1;			gb = p1->gb;					optimal = 0;					
	}

	//deconstructor
	~pattern(){
#ifdef SAVERESULTS
		if( stage == 1 || stage == 2 || stage == 3 ) delete M;
#else
		if( stage == 1 || stage == 2 || stage == 3 || stage == 4)       	delete M;
#endif
		
		//delete the coccs and occs
		delete coccs;
	}

	//member functions
	static void initPattern(gBase *gb1, MyTimer & t, Register & mr);
	void setChild(EDGES & cm, EDGES & nm);
	void fsm1( vector<AdjMatrix * > & result, int freq, int level, int density, 
		string nodef, string edgef, string header, int sizelimit, 
		int sizeup, map<int,int> &gIndex);
	void print();
	DLONG   getTotalPatterns()  { return totalPatterns; }

};

typedef vector<pattern *> CHLTYPE ;

#endif

